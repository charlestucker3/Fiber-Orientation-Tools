function [C, varargout] = lielens(Cf, Cm, vf, aspect, varargin)
%LIELENS     Lielens model for stifness of an aligned-fiber composite.
%      C = LIELENS(CF, CM, VF, ASPECT) uses the Lielens or double-inclusion model to
%      compute the stiffness tensor C of a composite with
%      unidirectional reinforcement. Cf and Cm are the stiffness tensors
%      for the fiber and matrix in 6x6 form, and C is returned in the same
%      form.  VF is the fiber volume fraction and ASPECT is the fiber
%      aspect ratio (length/diameter) or the ellipsoidal aspect ratio.
%
%      Cm and Cf must be isotropic.  C will be transversely isotropic
%      about the 1 axis. 
%
%  [C, BETA] = LIELENS(CF, CM, VF, ASPECT, ALPHAF, ALPHAM) also computes the
%      thermal stress tensor BETA in 6x1 contracted form.  ALPHAF and
%      ALPHAM are the thermal expansion coefficients of the fiber and
%      matrix.  ALPHAM must be scalar (i.e., the matrix has isotropic
%      thermal expansion).  If ALPHAF is scalar the fiber is isotropic.  If
%      ALPHAF is 6x1 then the fiber is anisotropic, in which case ALPHAF
%      must be transversely isotropic around the 1 axis.  The thermal
%      expansion tensor of the composite can be computed as
%      ALPHA = INV4(C) * R4 * BETA.

%      Reference: G. Lielens, P. Pirotte, A. Couniot, F. Dupret, and R.
%      Keunings. Prediction of thermo-mechanical properties for
%      compression-moulded composites. Composites A, 29:63–70, 1997.
 

I4 = diag([1,1,1,0.5,0.5,0.5]);  % 4th-order unit tensor
R4 = diag([1,1,1,2,  2,  2]);    % Used in contracted-notation products


% Fiber and matrix Poisson ratios
EngConM = C2eng(Cm);
Num     = EngConM(9);
EngConF = C2eng(Cf);
Nuf     = EngConF(9);

% Eshelby tensors for fiber and matrix
EshM = eshtens(aspect, Num);
EshF = eshtens(aspect, Nuf);

% Strain concentration tensors for Hashin-Shtrikman-Willis lower and upper
% bounds
Blower = inv4( I4 + EshM*R4*inv4(Cm)*R4*(Cf-Cm) ); % = Mori-Tanaka
Bupper =       I4 + EshF*R4*inv4(Cf)*R4*(Cm-Cf)  ;

% Interpolated strain concentration tensor
xi = (vf + vf^2)/2;
B = inv4( (1-xi)*inv4(Blower) + xi*inv4(Bupper) );

% Composite stiffness
C = ( (1-vf)*Cm + vf*Cf*R4*B ) * R4 * inv4( (1-vf)*I4 + vf*B );

if nargout >= 2
    if nargin >= 6  % Compute the thermal stress tensor of the composite
        % Thermal expansion tensor for the fiber
        if length(varargin{1}) == 1 % alphaf is a scalar
            alphaf = varargin{1} * [1,1,1,0,0,0]';
        elseif length(varargin{1}) == 6 % alphaf is a tensor
            alphaf = reshape(varargin{1}, 6, 1);
        else
            error('alphaf must be scalar or 6x1')
        end
        % Thermal stress tensors for the fiber and matrix
        betaf = Cf * R4 * alphaf;
        betam = Cm * R4 * (varargin{2} * [1,1,1,0,0,0]');
        % Thermal stress tensor for the composite, Eqn. (8.89)
        beta = (Cf - C) *R4* inv4(Cf - Cm) *R4* betam ...
             + (Cm - C) *R4* inv4(Cm - Cf) *R4* betaf;
        varargout{1} = beta;
    else
        error('Not enough input arguments to compute thermal expansion')
    end
end

return
