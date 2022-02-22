function [EngCon, varargout] = halpin (Ef, Nuf, Em, Num, vf, aspect, ...
                                       varargin)
%HALPIN      Halpin-Tsai equations.
%     [ENGCON]  = HALPIN(EF, NUF, EM, NUM, VF, ASPECT) calculates the
%     linear elastic engineering (technical) constants of a short-fiber
%     composite with aligned fibers using the Halpin-Tsai equations.
%
%     [ENGCON, C]  = HALPIN(EF, NUF, EM, NUM, VF, ASPECT) also returns the
%     stiffness tensor C in 6x6 matrix form.
%
%     [ENGCON, C, BETA]  = HALPIN(EF, NUF, EM, NUM, VF, ASPECT, ALPHAF,
%     ALPHAM) also returns the thermal stress tensor BETA in 6x1 form.
%     ALPHAF and ALPHAM give the thermal expansion of the fiber and matrix,
%     respectively, and must be present when calculating BETA. The thermal
%     expansion tensor of the composite can be calculated as
%        ALPHA = INV4(C) * R4 * BETA;
%     ALPHAM must be a scalar.  ALPHAF can be a scalar (isotropic fiber
%     properties) or a second-order tensor in 6x1 form.  If ALPHAF is a
%     tensor, it must be transversely isotropic about the 1 axis.  
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

if nargout >= 2
    % Convert engineering constants to C tensor
    C = eng2C(EngCon);
    varargout{1} = C; 
end

if nargout >= 3
    if nargin >= 8  % Compute the thermal stress tensor of the composite
        % Thermal expansion tensor for the fiber
        if length(varargin{1}) == 1 % alphaf is a scalar
            alphaf = varargin{1} * [1,1,1,0,0,0]';
        elseif length(varargin{1}) == 6 % alphaf is a tensor
            alphaf = reshape(varargin{1}, 6, 1);
        else
            error('alphaf must be scalar or 6x1')
        end
        % Stiffness tensors for the fiber and matrix
        Cf = iso2C(Ef, Nuf);
        Cm = iso2C(Em, Num);
        % Thermal stress tensors for the fiber and matrix
        R4 = diag([1,1,1,2,2,2]);
        betaf = Cf * R4 * alphaf;
        betam = Cm * R4 * (varargin{2} * [1,1,1,0,0,0]');
        % Thermal stress tensor for the composite, Eqn. (8.89)
        beta = (Cf - C) *R4* inv4(Cf - Cm) *R4* betam ...
             + (Cm - C) *R4* inv4(Cm - Cf) *R4* betaf;
        varargout{2} = beta;
    else
        error('Not enough input arguments to compute thermal expansion')
    end
end

return
