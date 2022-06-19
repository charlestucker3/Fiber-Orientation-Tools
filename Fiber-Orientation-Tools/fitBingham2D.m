function [T, psi, phi, A, A4] = fitBingham2D(A, npoints)
%[T, PSI, PHI] = FITBINGHAM2D(A, NPOINTS) finds the second-order tensor T
%    in 2-D that gives a Bingham distribution, psi(p) = exp(p*T*p) matching
%    the planar second-order orientation tensor A.  NPOINTS is the number
%    of points covering one half the unit circle used to integrate the
%    distribution function.  The vector PSI is the distribution function at
%    angles PHI.  A and T are 2x2 symmetric matrices.
%
%[T, PSI, PHI, A, A4] = fitBingham2D(A, NPOINTS) also returns the
%    second-order tensor A and the fourth-order tensor A4 corresponding to
%    the fitted Bingham distribution.  The output A may not exactly equal
%    the input A due to discretization errors.


Av  = tens2vec(A);       % Target A tensor in 3x1 vector form
Tv = [-1.84; -1.84; 0];  % Initial guess at the vector form of T

maxiter = 20;    % Maximum iterations
tol     = 1e-12;  % Solution tolerance
dTv     = 0.001; % Increment in T components for finite diff. derivatives

% Newton-Raphson iteration
for i = 1:maxiter
    [Avtrial] = BinghamA(Tv, npoints);  % A vector at this value of Tv
    f = Avtrial - Av;                   % We are seeking f = 0
    if norm(f) < tol
        % The solution has converged
        break
    end
    
    % Form matrix of derivatives, K = df/dTv
    K = zeros(3);
    for j = 1:3
        Tdiff    = Tv;
        Tdiff(j) = Tdiff(j) + dTv;
        [Avdiff] = BinghamA(Tdiff, npoints); % This function appears below
        K(:,j)   = (Avdiff - Avtrial) / dTv;
    end
    
    % Update the solution
    Tv = Tv - K\f;
end

if norm(f) < tol
%     fprintf('fitBingham2D converged in %i iterations\n', i)
    T = vec2tens(Tv);
    [Av, A4v, psi, phi] = BinghamA(Tv, npoints);
    A  = vec2tens(Av);
    A4 = vec2tens4(A4v);
else
    fprintf('NO CONVERGENCE in fitBingham2D\n')
    fprintf('norm(f) = %e after %i iterations\n', norm(f), i)
end

end

%% The function BinghamA is used only within this function
function [Av, varargout] = BinghamA(Tv, npoints)
%AV = BINGHAMA(TV, NPOINTS) returns the planar second-order tensor AV in 3x1
%    vector form for a Bingham distribution, PSI = exp(p*T*p).  TV is the
%    tensor T stored in 3x1 vector form.  NPOINTS is the number of
%    equally-spaced integration points in phi and controls the accuracy of
%    the integration used to find AV.
%
%[AV, A4V] = BINGHAMA(TV, NPOINTS) also returns the fourth-order tensor A4,
%    in 5x1 column-vector form.  See tens2vec4 for the arrangment of the
%    components, and use A4 = vec2tens4(A4V) to get the tensor in 3x3
%    matrix form.
%
%[AV, A4V, PSI, PHI] = BINGHAMA(TV, NPOINTS) returns NPOINTS distribution
%    function values PSI and the corresponding angles PHI.
%
%  See also: TENS2VEC4, VEC2TENS4.


phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);           % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;             % Delta(phi)

pp = [cos(phi).*cos(phi); ...  % Each column gives the components of the
      sin(phi).*sin(phi); ...  % pp tensor for one phi value.
      sin(phi).*cos(phi)]; 
  
psi = exp(Tv' * diag([1,1,2]) * pp); % Bingham distribution function

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