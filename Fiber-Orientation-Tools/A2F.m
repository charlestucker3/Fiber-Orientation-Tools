function [F, varargout] = A2F(A, F1func, F2func, varargin)
%F = A2F(A, F1func, F2func) finds the deformation tensor F (3x3) 
%    corresponding to the second-order orientation tensor A (3x3)
%    under the assumptions of the Natural closure approximation. 
%
%F = A2F(A, F1func, F2func, TOL) improves the solution using
%    Newton-Raphson iteration until the eigenvalues of the A recovered 
%    from F are within TOL of the input values.  If TOL is not given on
%    input, a value of 1e-3 is used.  This usually requires iterations
%    only very close to the U corner.  
%
%[F, NITER, AERR] = A2F(A, F1func, F2func, TOL) also returns
%    the number of Newton iterations NITER and AERR, the maximum absolute 
%    error in the eigenvalues of A.    
%
%    The Matlab scattered interpolant objects F1func and F2func 
%    must be passed to the function.  These are normally loaded in the 
%    calling program using LOAD A2Ffuncs.  

% --Find principal values and axes of A 
[Evals, R] = eigsort(A);
% [Evecs, Aprin] = eig(A);                 % eigenvectors/values of A
% [Evals, sortOrder] = sort(diag(Aprin));  % sort Evals in ascending order
% Evals = flipud(Evals);                   % change to decreasing order
% % Sorted eigenvectors form the rotation matrix
% R = [Evecs(:,sortOrder(3)), Evecs(:,sortOrder(2)), Evecs(:,sortOrder(1))];
%   Now A = R * diag(Evals) * R'

% -- Interpolate to get Fp, the F tensor in the principal axes of A
Fp1 = F1func(Evals(1), Evals(2));
Fp2 = F2func(Evals(1), Evals(2)); 
Fp3 = 1/(Fp1*Fp2); % This guarantees det(F) = 1

if isnan(Fp1) || isnan(Fp2)
    % This could occur if the query points are outside the UBT triangle
    fprintf('Fp1 and/or Fp2 are NaN.  Input A is:\n')
    fprintf('  %f   %f   %f\n', A')
    fprintf('Fp1 = %f,  Fp2 = %f\n', Fp1, Fp2)
    error('Cannot process NaN values of F')
end

%  - Set the error tolerance
if nargin >= 4
    tol = varargin{1};  % User-supplied value
else
    tol = 1e-3;   % The default tolerance
end

% -- Newton-Raphson iteration to improve the accuracy of the solution,
%    if desired (or needed).  Note that this is always an option, even
%    if the user has not supplied a tolerance, and some iterations are
%    often needed very close to the U corner.  
%    Check intial error
Af = F2A(diag([Fp1, Fp2, Fp3])); % A for this F, in principal axes
aerr = max(abs(Evals - diag(Af)));
%    Limit iterations, to be safe
niter   = 0;
maxiter = 20; % Can usually get errors ~1e-12 in 3 or 4 iterations,
%    but may need more very near the U corner.
%  - Main iteration loop
while aerr > tol && niter < maxiter
    % -- Find the Jacobian matrix by finite differences
    J = zeros(2,2);
    dfJ = 0.001;   % delta-F for Jacobian (somewhat arbitrary)
    Afd1 = F2A(diag([Fp1+dfJ, Fp2,    1./((Fp1+dfJ)*Fp2)]));
    Afd2 = F2A(diag([Fp1   , Fp2+dfJ, 1./(Fp1*(Fp2+dfJ))]));
    J(1,1) = (Afd1(1,1) - Af(1,1)) / dfJ;
    J(2,1) = (Afd1(2,2) - Af(2,2)) / dfJ;
    J(1,2) = (Afd2(1,1) - Af(1,1)) / dfJ;
    J(2,2) = (Afd2(2,2) - Af(2,2)) / dfJ;
    % -- Update the F values
    dF = J \ [Evals(1)-Af(1,1); Evals(2)-Af(2,2)];
    Fp1 = Fp1+dF(1); Fp2 = Fp2+dF(2); Fp3 = 1/(Fp1*Fp2);
    % -- Get A for the updated F, and computer the error
    Af = F2A(diag([Fp1, Fp2, Fp3]));
    aerr = max(abs(Evals - diag(Af)));
    niter = niter + 1;
end  % Normally will iterate to convergence
%     fprintf('A2Fnewton: aerr = %8.3e in %i iterations\n', aerr, niter)
if niter == maxiter && aerr > tol
    fprintf('WARNING: A2F did not converge in %i iterations\n', niter)
    fprintf('         aerr = %8.3e\n', aerr)
end


% -- Return F to the original coordinate system
F = R * diag([Fp1, Fp2, Fp3]) * R';

% -- Fill in any variable output arguments
if nargout >= 2
    varargout{1} = niter;
end
if nargout >= 3
    varargout{2} = aerr;
end



return
