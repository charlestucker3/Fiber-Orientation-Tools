function [T, psi, phi, A, A4] = fitERdistn2D(A, npoints)
%[T, PSI, PHI, A, A4] = FITERDISTN2D(A, NPOINTS) finds the 
%    second-order tensor T in 2-D that gives a ellipse-radius distribution,
%    psi(p) = (p*T*p)^(-1/2r) matching the planar second-order orientation
%    tensor A. NPOINTS is the number of points covering one half the unit
%    circle used to integrate the distribution function.  The vector PSI is
%    the distribution function at angles PHI.  A and T are 2x2 symmetric
%    matrices.
%
%    See also: ERdistnA

% Eigenvalues and eigenvectors of A
[lambda, R] = eigsort(A);
if lambda(1) > 0.9999
    error('fitERdistn2D: lambda1 = %8.6f is greater than 0.9999', lambda(1))
end

% Find the ER parameter r in the principal axis system of A.
% The function r2lambda1 is defined at the bottom of this file.
r = fzero(@(r) r2lambda1(r, npoints) - lambda(1), [0.007, 1.1]);

% Build the T tensor
T = diag([r^2, 1]);              % T in principal axes, not normalized
phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);             % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;               % Delta(phi)
% psi  = ((1-r^2) * sin(phi).^2 + r^2).^(-1/(2*r)); % ODF kernel
psi  = (r^2 * cos(phi).^2 + sin(phi).^2).^(-1/(2*r)); % ODF kernel
Z    = 1 / (2*sum(psi) * dphi) ; % Normalization factor
T = Z^(-2*r) * T;                % Normalized T
T = R * T * R';                  % Return T to laboratory axes

% Compute additional information
% The function ERdistnA appears at the bottom of the file
[Av, A4v, psi, phi] = ERdistnA(T, npoints);
A  = vec2tens(Av);
A4 = vec2tens4(A4v);

end

%% ------------------------------------------------------------

function lambda1 = r2lambda1(r, npoints)
%LAMBDA1 = R@LAMBDA1(R, NPOINTS) finds the largest eigenvalue of the
%    second-order orientation tensor for the plnar ellipse-radius (ER)
%    distribution with parameter R.  NPOINTS is the number of equally
%    spaced points on half a circle used for the integration.

phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);           % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;             % Delta(phi)

psi     = ((1-r^2) * sin(phi).^2 + r^2).^(-1/(2*r)); % ODF kernel
psi     = psi / (2*sum(psi) * dphi) ;                % Normalized ODF
lambda1 = 2*sum(psi .* cos(phi).^2) * dphi;          % <cos^2 phi>

end

%% ------------------------------------------------------------

function [Av, varargout] = ERdistnA(T, npoints)
%AV = ERDISTNA(T, NPOINTS) returns the planar second-order tensor AV 
%    in 3x1 vector form for a elllipse-radius (ER) distribution, PSI =
%    (p*T*p)^(-1/2r).  T is a tensor stored in 2x2 form. NPOINTS is the
%    number of equally-spaced integration points in phi, and controls the
%    accuracy of the integration used to find AV.
%
%[AV, A4V] = ERDISTNA(T, NPOINTS) also returns the fourth-order 
%    tensor A4V in 5x1 column-vector form.  See tens2vec4 for the
%    arrangment of the components.  Use A4 = vec2tens4(A4V) to get the
%    fourth-order tensor in 3x3 matrix form.
%
%[AV, A4V, PSI, PHI] = ERDISTNA(T, NPOINTS) also returns NPOINTS 
%    distribution function values PSI and the corresponding angles PHI.
%
%  See also: FITERDISTN2D, TENS2VEC4, VEC2TENS4.


phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);           % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;             % Delta(phi)

pp = [cos(phi).*cos(phi); ...  % Each column gives the components of the
      sin(phi).*sin(phi); ...  % pp tensor for one phi value.
      sin(phi).*cos(phi)]; 
  
% r is the ratio b/a of the minor and major axis lengths of the ellipse
[Eval, ~] = eigsort(T);     % Eigenvalues of T, sorted
r = sqrt(Eval(2)/Eval(1));  % r parameter for the ER distribution
  
psi = (tens2vec(T)' * diag([1,1,2]) * pp).^(-1/(2*r)); % Ellipse-radius dist'n.

Av = 2 * pp * psi' * dphi;     % Second-order tensor, contracted form

if nargout >= 2
    % Calculate the fourth-order tensor
    p4 = [cos(phi).^4; ...     % Each column of p4 contains the components
          sin(phi).^4; ...     % of the pppp tensor for one phi value.
          cos(phi).^2 .* sin(phi).^2; ...
          cos(phi).^3 .* sin(phi); ...
          cos(phi)    .* sin(phi).^3];

      A4v = 2 * p4 * psi' * dphi;     % Fourth-order tensor, contracted
      varargout{1} = A4v;
end

% Return additional arguments if requested
if nargout >= 3
    varargout{2} = psi;
end
if nargout >= 4
    varargout{3} = phi;
end


end