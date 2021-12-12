function tau = tauFiber(A, D, etaI, Np, closure)
%TAU = TAUFIBER(A, D, ETAI, NP, CLOSURE) finds the extra-stress tensor TAU
%      for a fiber suspension with second-order orientation tensor A (3x3),
%      rate-of-deformation tensor D (3x3), isotropic viscosity ETAI and
%      particle number NP, using the closure approximation specified by
%      CLOSURE.
%
%      Options for CLOSURE are as in closeA4.m (e.g., 'N' for natural, 'H'
%      for hybrid, 'I' for IBOF).  

Afour = closeA4(A, closure);             % Fourth-order tensor
R4 = diag([1, 1, 1, 2, 2, 2]);           % Used in A4:D
Dv = tens2vec(D);                        % D in 6x1 form
tau = 2 * etaI * (Dv + Np*Afour*R4*Dv);  % The constitutive eqn. 

return