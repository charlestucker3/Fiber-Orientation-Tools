function [diffParam, Av, varargout] = fitARD(A11, A33, xi, closure, ...
                                             diffModel, diffParam, varargin)
%[DIFFPARAM, AV] = FITARD(A11, A33, XI, CLOSURE, DIFFMODEL, DIFFPARAM) 
%      finds the fiber orientation model parameters DIFFPARAM for the 
%      ARD model DIFFMODEL that match the orientation components A11 and 
%      A33  in 1-3 steady simple shear flow, v(1) = G * x(3).  The particle
%      shape factor is XI and CLOSURE gives the closure approximation that
%      is used.  The input parameter DIFFPARAM provides an initial guess
%      for the model parameters.  AV is the steady-state orientation, in
%      6x1 column-vector form.
%
%[DIFFPARAM, AV] = FITARD(A11, A33, XI, CLOSURE, DIFFMODEL, ...
%      DIFFPARAM, KINMODEL, KINPARAM) uses the orientation-kinetics model
%      KINMODEL with parameters KINPARAM.  If these are absent, the default
%      is Jeffery kinetics.  See Adot2.m for the available values of
%      KINMODEL and KINPARAM.
%
%[DIFFPARAM, AV, STABLE] = FITARD(... also returns STABLE = 1 if the
%      steady-state solution is stable, and STABLE = 0 if the steady state
%      is unstable.  This usually occurs when the requested (A11, A33) pair
%      is outside the valid range of the ARD model.  STABLE = -1 if the
%      solution did not converge.  
%
%      The function uses Newton-Raphson to find dA/dt = 0 for A11, A33
%      and A13.  A13 is adjusted, along with two model parameters, to
%      achieve this.  If the MRD model is used, D2 is held fixed at its
%      input value.  
% 
%      CLOSURE, DIFFMODEL and DIFFPARAM have the same meanings as in
%      ADOT2.m.  In addition, for the Wang model DIFFPARAM(3:4) give the 
%      alternate parameters, CI and w.

L = zeros(3,3); L(1,3) = 1;  % Velocity gradient, 1-3 simple shear

% Shell orientation to be matched, 6x1 vector form
Av = zeros(6,1);
Av(1) = A11;
Av(3) = A33;
Av(2) = 1 - A11 - A33;
Av(5) = 0.10;  % Arbitrary initial value; adjusted during sol'n.

% Solution tolerances
maxiter = 20;    % Max iterations for Newton-Raphson
errTol = 1e-14;  % Tolerance on dAdt = 0.

% Orientation kinetics model and parameters
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
% The x vector contains the independent variables we are solving for
switch upper(diffModel)    
    case {'W', 'I', 'P'} 
        % I = iARD model:        x = [CI, CM,    A13]';
        % P = pARD model:        x = [CI, Omega, A13]';
        % W = Wang (WPT) model:  x = [b1, b3,    A13]';
         x = [diffParam(1), diffParam(2), Av(5)]';
         % Extend diffParam for compatibility in calling Adot2.m
         diffParam = [diffParam(1:2), 0];  
    case 'M'  
        % M = MRD model:  x = [CI, D3, A13]';  D2 = input value
         x = [diffParam(1), diffParam(3), Av(5)]';
    otherwise
        error('diffModel = %s is not supported', diffModel)
end

% The y vector contains the dependent variables whose value we want = 0
dAdt = Adot2(Av, L, xi, closure, diffModel, diffParam, kinModel, kinParam);
y = dAdt([1,3,5]);
% % yerr = sum(abs(y));  % Scalar measure of convergence
yerr = norm(y);


%%  - Main iteration loop
niter   = 0;
while yerr > errTol && niter < maxiter
    % -- Find the Jacobian matrix by finite differences
    J = zeros(3,3);
    delta = 0.001;   % delta for model params or A13 (somewhat arbitrary)
    dAdt1 = Adot2(Av, L, xi, closure, diffModel, diffParam+[delta,0,0], ...
                  kinModel, kinParam);
    ydx1 = dAdt1([1,3,5]);
    J(:,1) = (ydx1 - y) / delta;
    switch upper(diffModel)
        case {'W', 'P', 'I'}
    dAdt2 = Adot2(Av, L, xi, closure, diffModel, diffParam+[0,delta,0],...
                  kinModel, kinParam);
        case 'M'
    dAdt2 = Adot2(Av, L, xi, closure, diffModel, diffParam+[0,0,delta],...
                  kinModel, kinParam);
        otherwise
            error('Unexpected value of diffModel = %s', diffModel)
    end
    ydx2 = dAdt2([1,3,5]);
    J(:,2) = (ydx2 - y) / delta;
    dAdt3 = Adot2(Av + [0,0,0,0,delta,0]', L, xi, closure, ...
                  diffModel, diffParam, kinModel, kinParam);
    ydx3 = dAdt3([1,3,5]);
    J(:,3) = (ydx3 - y) / delta;
    % -- Update the x values using Newton-Raphson 
    dx = -J\y;
    x = x + dx;
    % --  Update diffParam and Av
    switch upper(diffModel)
        case {'W', 'P', 'I'}
            diffParam = [x(1), x(2), 0];
        case 'M'  % Note that D2 from the input is re-used
            diffParam = [x(1), diffParam(2), x(2)];
        otherwise
            error('Unexpected value of diffModel = %s', diffModel)
    end
    Av(5) = x(3);  % Update A13
    % -- Update y and compute the overall error
    dAdt = Adot2(Av, L, xi, closure, diffModel, diffParam, ...
                 kinModel, kinParam);
    y = dAdt([1,3,5]);
% %     yerr = sum(abs(y));
    yerr = norm(y);
    niter = niter + 1;
end  % Normally will iterate to convergence

if niter == maxiter && yerr > errTol
    steady = -1;
    fprintf('WARNING: fitARD did not converge in %i iterations\n', niter)
    fprintf('         yerr = %8.3e\n', yerr)
else
    % Check the stability of the solution.  Use Asteady to compute the
    % Jacobian for d(Adot)/dA and make sure all eigenvalues have negative
    % real parts.  Asteady.m will print a warning if Av is an unstable
    % steady state.
    [~, steady] = Asteady(Av, L, xi, closure, diffModel, diffParam, ...
        kinModel, kinParam);
end
if nargout >= 3
    varargout{1} = steady;
end

% Return diffParam to its final size for the two-parameter models
switch upper(diffModel)
    case {'I','P'}
        diffParam = diffParam(1:2);
    case {'W'}
        diffParam = [diffParam(1:2), 0, 0];  % Add two entries
        diffParam(3) = diffParam(1) + diffParam(2);  % CI for Wang
        diffParam(4) = diffParam(2) / diffParam(3);  % w for Wang
    case {'M'}
       % diffParam is already the right size
end

% Check parameter values and issue a warning if any are invalid.
valid = true(1);  % Flag for invalid model parameters
if diffParam(1) <= 0 % CI or b1 should be >= 0 for all models
    valid = false(1);
end
switch upper(diffModel)
    case 'I'
        if diffParam(2) < 0 || diffParam(2) > 1  % Want 0 <= CM <= 1
            valid = false(1);
        end            
    case 'P'
        if diffParam(2) < 0.5 || diffParam(2) > 1 % Want 0 <= Omega <= 0.5
            valid = false(1);
        end
    case 'M'
        if diffParam(3) < 0   % Want D3 >= 0
            valid = false(1);
        end
    case 'W'
        if diffParam(2) < 0 || ...   % Want 0 <= w <= 1
           diffParam(2)/(diffParam(1) + diffParam(2)) > 1 
            valid = false(1);
        end
end
if ~valid
  fprintf('WARNING: fitARD found illegal model parameters\n')
end



return
