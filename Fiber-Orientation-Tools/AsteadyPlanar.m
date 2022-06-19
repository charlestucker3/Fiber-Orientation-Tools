function [Av, varargout]  = AsteadyPlanar(Av, L, xi, closure, ...
                                          diffModel, diffParam)
%AV = ASTEADYPLANAR(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM) finds the 
%      planar steady-state orientation tensor Av (2x1 vector form)
%      corresponding to the velocity gradient L (2x2), particle shape
%      factor XI, closure approximation CLOSURE, and rotary diffusion model
%      DIFFMODEL with parameter DIFFPARAM.  These input parameters have the
%      same meaning as in AdotPlanar.m.  AV on input is an initial guess at
%      the steady orientation state.  The method is Newton-Raphson, so a
%      reasonable initial guess may be needed for the function to converge.
%
%      Note that slow kinetics models like SRF and RSC do not affect the
%      steady-state solution, so they are not needed as input.
%
%[AV, STABLE] = ASTEADYPLANAR(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM) 
%      also returns information about the stability of the solution. STABLE
%      = 1 if the solution is stable (all eigenvalues of d(Adot)/dA have
%      negative real parts) and 0 if not.  STABLE = -1 if the solution did
%      not coverge.
%
%      See also: ADOTPLANAR, CLOSEA4PLANAR

% Ensure that Av is a column vector
if isrow(Av)
    Av = Av';
end

% Set the rotary diffusion model and interaction coefficient
if strcmpi(diffModel, 'J')
    CI = 0;   % No interaction coefficient for the Jeffery model
              % **** WARNING: This many not have a steady state! ****
elseif strcmpi(diffModel, 'F')
    % Folgar-Tucker rotary diffusion
    CI = diffParam;
else
    error('%s is not a legal value for diffModel', diffModel)
end

% Solution tolerances
maxiter = 50;    % Max iterations for Newton-Raphson
errTol = 1e-12;   % Tolerance on dAdt = 0.

% Scale L so that dA/dt values will be well scaled.
L = L / max(max(abs(L)));
                    
                    
%% Initial guess at independent variables, and initial trial solution
% The x vector contains the independent variables we are solving for, which
% are A11 and A13.  The A22 component is found from A11 by normalization.

x = Av([1,3]);
% The y vector contains the dependent variables whose value we want = 0
dAdt = AdotPlanar(Av, L, xi, closure, diffModel, CI);
y = dAdt([1,3]);
yerr = sum(abs(y));  % Scalar measure of convergence
yerrplot = zeros(maxiter, 1);
yerrplot(1) = yerr;


%%  - Main iteration loop
niter   = 0;
while yerr > errTol && niter < maxiter
    % -- Form the Jacobian matrix by finite differences
    J = zeros(2,2);
    delta = 0.001;   % delta for tensor components (somewhat arbitrary)
    dAdt = AdotPlanar(Av + [delta,-delta,0]', L, xi, closure, ...
                      diffModel, CI);
    J(:,1) = (dAdt([1,3]) - y) / delta;
    dAdt = AdotPlanar(Av + [0,0,delta]',      L, xi, closure, ...
                     diffModel, CI);
    J(:,2) = (dAdt([1,3]) - y) / delta;
    
    % -- Update the x vector using Newton-Raphson 
    dx = -J\y;
    x = x + dx;
    % --  Update Av
    Av = [x(1), 1-x(1), x(2)]';
    % -- Update y and compute the new overall error
    dAdt = AdotPlanar(Av, L, xi, closure, diffModel, CI);
    y = dAdt([1,3]);
    yerr = sum(abs(y));
    niter = niter + 1;
    yerrplot(niter+1) = yerr;
end  % Normally will iterate to convergence


if niter == maxiter && yerr > errTol
    stable = -1;
    fprintf('WARNING: AsteadyPlanar did not converge in %i iterations\n', niter)
    fprintf('         yerr = %8.2e  diffParam = %8.2e\n', yerr, diffParam)
else
    % Check the stability if the solution.  This requires re-computing J,
    % since the iteration loop may have been skipped if the initial guess
    % was also the final answer.
    J = zeros(2,2);
    delta = 0.001;   % delta for tensor components (somewhat arbitrary)
    dAdt = AdotPlanar(Av + [delta,-delta,0]', L, xi, closure, ...
                      diffModel, CI);
    J(:,1) = (dAdt([1,3]) - y) / delta;
    dAdt = AdotPlanar(Av + [0,0,delta]',      L, xi, closure, ...
                     diffModel, CI);
    J(:,2) = (dAdt([1,3]) - y) / delta;
    
    if max(real(eig(J))) >= 0
        fprintf('WARNING: AsteadyPlanar found an unstable steady state\n')
        fprintf('    CI = %8.2e, closure = %s\n', CI, closure)
        fprintf('    L = [%5.2f, %5.2f; %5.2f, %5.2f]\n', L')
        stable = 0;
    else
        stable = 1;
    end
end

if nargout >= 2
    varargout{1} = stable;
end

return
