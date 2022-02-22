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

% Use fzero and A11steady (see below) to find CI
CI = fzero(@(Ctrial) A11 - A11steady(L, xi, closure, 'F', Ctrial), [1e-6 1]);

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


%% ---------------------------------------------------------- %%

function A11 = A11steady(L, xi, closure, diffModel, diffParam)
% A11 steady is used only by fitCI.   It calls Asteady, and possibly
% Adot2, to find the steady-state value of A11 in 6x1 vector form.

Aviso = tens2vec(eye(3)/3);  % Isotropic orientation, vector form

% Set initial guess for Asteady
if (closure == 'I' || closure == 'H') && diffParam < 1e-3
    % Need to integrate the ODE for a while to get close to steady state
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    tspan = 0:20;
    [~, Av] = ode45(@(t,Av) Adot2(Av, L, xi, closure, diffModel, diffParam), ...
                    tspan, Aviso, options);
    Avzero = Av(end,:)';
else
    Avzero = Aviso; % The isotropic state will work for most cases
end

% Find the steady-state A tensor
Av = Asteady(Avzero, L, xi, closure, diffModel, diffParam);
A11 = Av(1);

return
end
