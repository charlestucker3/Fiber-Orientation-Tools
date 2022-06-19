function [T4, psi, phi, A, A4] = fitMaxEntropy2D(A4, npoints)
%[T4, PSI, PHI, A, A4] = FITMAXENTROPY2D(A4, NPOINTS) finds the 
%    fourth-order tensor T4 in 2-D that gives a fourth-order
%    maximum-entropy distribution, psi(p) = exp(pp:T4:pp) matching the
%    planar fourth-order orientation tensor A4.  NPOINTS is the number of
%    points covering one half the unit circle used to integrate the
%    distribution function.  The vector PSI is the distribution function at
%    angles PHI, and A and A4 are the orientation tensors for the fitted
%    function.  A4 and T4 are stored in contracted form as 3x3 symmetric
%    matrices.
%
%    The second-order maximum entropy distribution is found using
%    fitBingham2D.
%
%    See also: FITBINGHAM2D

A4v  = tens2vec4(A4);   % Target A4 tensor in 5x1 vector form

T4v = [-1.84; -1.84; -0.61; 0; 0];  % Initial guess at the vector form of Tv

maxiter = 100;    % Maximum iterations
tol     = 1e-12;  % Solution tolerance
dTv     = 0.0001; % Increment in T4 components for finite diff. derivatives

% Newton-Raphson iteration
for i = 1:maxiter
    [~, A4vtrial] = maxEntropyA(T4v, npoints);  % A4v at this value of T4v
    f = A4vtrial - A4v;                      % We are seeking f = 0
    if norm(f) < tol
        % The solution has converged
        break
    end
    
    % Form matrix of derivatives, K = df/dT4v
    K = zeros(5);
    for j = 1:5
        Tdiff    = T4v;
        Tdiff(j) = Tdiff(j) + dTv;
        % The function maxEntropyA appears at the bottom of this file
        [~, A4vdiff]   = maxEntropyA(Tdiff, npoints); 
        K(:,j)   = (A4vdiff - A4vtrial) / dTv;
    end
    
    % Update the solution 
    T4v = T4v - K\f;
end

if norm(f) < tol
%     fprintf('fitMaxEntropy2D converged in %i iterations\n', i)
else
    fprintf('NO CONVERGENCE. norm(f) = %e after %i iterations\n', norm(f), i)
end

% Return current results, even with no convergence.
T4 = vec2tens4(T4v);
[Av, A4v, psi, phi] = maxEntropyA(T4v, npoints);
A  = vec2tens(Av);
A4 = vec2tens4(A4v);

end


%% The function maxEntropyA is used only within this program

function [Av, varargout] = maxEntropyA(T4v, npoints)
%AV = MAXENTROPYA(T4V, NPOINTS) returns the planar second-order tensor AV 
%    in 3x1 vector form for a fourth-order maximum entropy distribution,
%    PSI = exp(pp:T4:pp).  T4V is the tensor T4 stored in 5x1 vector form.
%    NPOINTS is the number of equally-spaced integration points in phi and
%    controls the accuracy of the integration used to find AV.
%
%[A4V, AV] = MAXENTROPYA(T, NPOINTS) also returns the fourth-order tensor 
%    A4V, in 5x1 column-vector form.  See tens2vec4 for the arrangment of
%    the components, and use A4 = vec2tens4(A4V) to get the tensor in 3x3
%    matrix form. 
%
%[AV, A4V, PSI, PHI] = MAXENTROPYA(TV, NPOINTS) returns NPOINTS 
%    distribution function values PSI and the corresponding angles PHI.
%
%  See also: TENS2VEC4, VEC2TENS4.


phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);           % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;             % Delta(phi)

p4 = [cos(phi).^4; ...     % Each column of p4 contains the components
      sin(phi).^4; ...     % of the pppp tensor for one phi value.
      cos(phi).^2 .* sin(phi).^2; ...
      cos(phi).^3 .* sin(phi); ...
      cos(phi)    .* sin(phi).^3];
  
psi = exp(T4v' * diag([1,1,6,4,4]) * p4); % Max entropy dist'n. function


pp = [cos(phi).*cos(phi); ...  % Each column gives the components of the
      sin(phi).*sin(phi); ...  % pp tensor for one phi value.
      sin(phi).*cos(phi)]; 

Av = 2 * pp * psi' * dphi;     % Second-order tensor, contracted form

if nargout >= 2
    % Return the fourth-order tensor
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

