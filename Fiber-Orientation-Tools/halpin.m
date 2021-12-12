function [EngCon, varargout] = halpin (Ef, Nuf, Em, Num, vf, aspect)
%HALPIN      Halpin-Tsai equations.
%     [ENGCON]  = HALPIN(EF, NUF, EM, NUM, VF, ASPECT) calculates the
%     linear elastic eingineering (technical) constants of a short-fiber
%     composite with aligned fibers using the Halpin-Tsai equations.
%
%     [ENGCON, C]  = HALPIN(EF, NUF, EM, NUM, VF, ASPECT) also returns the
%     stiffness tensor C in 6x6 matrix form.
%
%     Input:
%          EF     = Fiber Young's modulus
%          NUF    = Fiber Poisson ratio
%          EN     = Matrix Young's modulus
%          NUM    = Matrix Poisson ratio
%          VF     = Fiber volume fraction
%          ASPECT = Fiber aspect ratio (L/D)
%
%     Output: 
%          ENGCON(1) = E11
%          ENGCON(2) = E22
%          ENGCON(3) = E33
%          ENGCON(4) = G23
%          ENGCON(5) = G31
%          ENGCON(6) = G12
%          ENGCON(7) = Nu23
%          ENGCON(8) = Nu31
%          ENGCON(9) = Nu12
%

% Storage for the results
EngCon = zeros(9,1);

% Modulus in the fiber direction, E11
zeta = 2*aspect;
eta = ((Ef/Em)-1) / ((Ef/Em)+zeta);
EngCon(1) = Em * (1 + zeta*eta*vf) / (1 - eta*vf);

% Major Poisson ratio, Nu12
EngCon(9) = Nuf*vf + Num*(1-vf);

% Modulus in the transverse direction, E22 and E33
zeta = 2;
eta = ((Ef/Em)-1) / ((Ef/Em)+zeta);
EngCon(2) = Em * (1 + zeta*eta*vf) / (1 - eta*vf);
EngCon(3) = EngCon(2);

% Axial shear modulus, G12 and G31
Gf = Ef / (2*(1+Nuf));
Gm = Em / (2*(1+Num));
zeta = 1;
eta = ((Gf/Gm)-1) / ((Gf/Gm)+zeta);
EngCon(6) = Gm * (1 + zeta*eta*vf) / (1 - eta*vf);
EngCon(5) = EngCon(6);

% Transverse shear modulus, G23
% Uses continuous-fiber version
zeta = (1+Num) / (3 - Num - 4*Num^2);
eta = ((Gf/Gm)-1) / ((Gf/Gm)+zeta);
EngCon(4) = Gm * (1 + zeta*eta*vf) / (1 - eta*vf);

% -- Dependent properties:
%    Transverse Poisson ratio, Nu23
EngCon(7) = EngCon(2)/(2.0*EngCon(4)) - 1.0;
%    Poisson ratio Nu31:
EngCon(8) = EngCon(9) * EngCon(3) / EngCon(1);

if nargout > 1
    % Convert engineering constants to C tensor
    varargout{1} = eng2C(EngCon);
end
%
return
