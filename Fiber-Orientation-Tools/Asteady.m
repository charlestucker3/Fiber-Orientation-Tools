function [Av, varargout]  = Asteady(Av, L, xi, closure, ...
                                    diffModel, diffParam, varargin)
%AV = ASTEADY(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM) finds the
%      steady-state orientation tensor Av (6x1 vector form) corresponding
%      to the velocity gradient L (3x1), particle shape factor XI, closure
%      approximation CLOSURE, and diffusion model DIFFMODEL with
%      parameter(s) DIFFPARAM.  AV on input is an initial guess at the
%      steady orientation state.  The method is Newton-Raphson, so a
%      reasonable initial guess may be needed for the function to converge.
%
%AV = ASTEADY(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM, KINMODEL, KINPARAM)
%      uses the orientation kinetics model KINMODEL with parameter
%      KINPARAM.  If these arguments are absent, Jeffery kinetics are used.
%      The steady-state solution is actually independent of the kinetic
%      model and its parameters, so this option mainly serves to confirm
%      that numerically.  
%
%[AV, STABLE] = ASTEADY(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM) also
%      returns information about the stability of the solution. STABLE = 1
%      if the solution is stable (all eigenvalues of d(Adot)/dA have
%      negative real parts) and 0 if not.  STABLE = -1 if the solution did
%      not coverge.
%
%      CLOSURE, DIFFMODEL, DIFFPARAM, KINMODEL and KINPARAM have the same
%      meanings as in Adot2.m


% Solution tolerances
maxiter = 150;    % Max iterations for Newton-Raphson
errTol = 1e-12;  % Tolerance on dAdt = 0.

% Scale L so that dA/dt values will be well scaled.
L = L / max(max(abs(L)));

% Set orientation kinetics model and parameters
if nargin >= 8
    % Use the values from the input arguments
    kinModel = varargin{1};
    kinParam = varargin{2};
else
    % Use default values
    kinModel = 'J'; % Jeffery (standard) kinetics
    kinParam =  1;
end
                    
                    
%% Initial guess at independent variables, and initial trial solution
% The x vector contains the independent variables we are solving for,
% which are A11, A33, A23, A31, A12.  The A22 component is found from A11
% and A33 by normalization.

x = Av([1,3:6]);
% The y vector contains the dependent variables whose value we want = 0
dAdt = Adot2(Av, L, xi, closure, diffModel, diffParam, kinModel, kinParam);
y = dAdt([1,3:6]);
yerr = sum(abs(y));  % Scalar measure of convergence


%%  - Main iteration loop
niter   = 0;
while yerr > errTol && niter < maxiter
    % -- Form the Jacobian matrix by finite differences
    J = zeros(5,5);
    delta = 0.001;   % delta for tensor components (somewhat arbitrary)
    dAdt = Adot2(Av + [delta,-delta,0,0,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,1) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,-delta,delta,0,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,2) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,delta,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,3) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,0,delta,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,4) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,0,0,delta]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,5) = (dAdt([1,3:6]) - y) / delta;
    
    % -- Update the x vector using Newton-Raphson 
    dx = -J\y;
    % **** Try some relaxation, to help with tricky cases ***
    if niter <= 50
        x = x + 0.1*dx;
    else
    x = x + dx;
    end
    % ***********************
    % --  Update Av
    Av = [x(1), 1-x(1)-x(2), x(2), x(3), x(4), x(5)]';
    % -- Update y and compute the new overall error
    dAdt = Adot2(Av, L, xi, closure, diffModel, diffParam, ...
                 kinModel, kinParam);
    y = dAdt([1,3:6]);
    yerr = sum(abs(y));
    niter = niter + 1;
end  % Normally will iterate to convergence

if niter == maxiter && yerr > errTol
    stable = -1;
    fprintf('WARNING: Asteady did not converge in %i iterations\n', niter)
    fprintf('         yerr = %8.2e  diffParam = %8.2e\n', yerr, diffParam)
else
    % Check the stability if the solution.  This requires re-computing J,
    % since the iteration loop may have been skipped if the initial guess
    % was also the final answer.
    J = zeros(5,5);
    delta = 0.001;   % delta for tensor components (somewhat arbitrary)
    dAdt = Adot2(Av + [delta,-delta,0,0,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,1) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,-delta,delta,0,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,2) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,delta,0,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,3) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,0,delta,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,4) = (dAdt([1,3:6]) - y) / delta;
    dAdt = Adot2(Av + [0,0,0,0,0,delta]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    J(:,5) = (dAdt([1,3:6]) - y) / delta;
    
    if max(real(eig(J))) >= 0
        fprintf('WARNING: Asteady found an unstable steady state\n')
        stable = 0;
    else
        stable = 1;
    end
end

if nargout >= 2
    varargout{1} = stable;
end

return
