function [CI, varargout] = fitCI(A11, xi, closure)
%CI = FITCI(A11, XI, CLOSURE) finds the value of the interaction
%    coefficient CI for the Folgar-Tucker model that gives a flow-direction
%    orientation tensor component of A11 at steady state in simple shear
%    flow.  XI is the particle shape factor (usually XI = 1) and CLOSURE is
%    the closure approximation, a single character as in closeA4.  
%
%[CI, A33, A13] = FITCI(A11, XI, CLOSURE) also returns the tensor
%    components A33 and A13 at steady state.  
%    The flow is v_1 = gammaDot * x_3.

% Keep target A11 values reasonable to avoid numerical failures
if A11 > 0.99
    error('A11 must be <= 0.99')
end
if A11 < 0.34
    error('A11 must be > 0.34')
end

% Velocity gradient: 1-3 simple shear
L = zeros(3); L(1,3) = 1;

% Use fzero to find CI
CI = fzero(@(Ctrial) A11-A11steady(L, xi, closure, 'F', Ctrial), [1e-6 1]);

% Find A13 and A33 if desired
if nargin >= 2
    Av = Asteady(tens2vec(eye(3)/3), L, xi, closure, 'F', CI);
    varargout{1} = Av(3);     % A33
    if nargin >= 3
        varargout{2} = Av(5); % A13
    end
end

return
end
