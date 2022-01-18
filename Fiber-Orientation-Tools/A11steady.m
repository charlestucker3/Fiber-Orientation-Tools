function A11 = A11steady(L, xi, closure, diffModel, diffParam)
% A wrapper for Asteady.m which returns only the A11 component.
% Used by fitCI.m.  See Asteady.m for an explanation of the arguments.

% Steady-state orientation, 6x1 vector form
Av = Asteady(tens2vec(eye(3)/3), L, xi, closure, diffModel, diffParam);
A11 = Av(1);

return
end