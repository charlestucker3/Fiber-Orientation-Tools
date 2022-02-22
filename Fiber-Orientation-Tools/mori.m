function [C, varargout] = mori(Cf, Cm, vf, aspect, varargin)
%MORI     Mori-Tanaka model for stifness of an aligned-fiber composite.
%   C = MORI(CF, CM, VF, ASPECT) uses the Eshelby/Mori-Tanaka model to
%      compute the stiffness tensor C of a composite with unidirectional
%      reinforcement. Cf and Cm are the stiffness tensors for the fiber and
%      matrix in 6x6 form, and C is returned in the same form.  VF is the
%      fiber volume fraction and ASPECT is the fiber aspect ratio
%      (length/diameter) or the ellipsoidal aspect ratio.
%
%      Cm must be isotropic, and Cf can be transversely isotropic about
%      the 1 axis.  C will be transversely isotropic about the 1 axis.
%
%  [C, BETA] = MORI(CF, CM, VF, ASPECT, ALPHAF, ALPHAM) also computes the
%      thermal stress tensor BETAn 6x1 contracted form.  ALPHAF and
%      ALPHAM are the thermal expansion coefficients of the fiber and
%      matrix.  ALPHAM must be scalar (i.e., the matrix has isotropic
%      thermal expansion).  If ALPHAF is scalar the fiber is isotropic.  If
%      ALPHAF is 6x1 then the fiber is anisotropic, in which case ALPHAF
%      must be transversely isotropic around the 1 axis.  The thermal
%      expansion tensor of the composite can be computed as
%      ALPHA = INV4(C) * R4 * BETA.
 

I4 = diag([1,1,1,0.5,0.5,0.5]);  % 4th-order unit tensor
R4 = diag([1,1,1,2,  2,  2]);    % Used in contracted-notation products


% Matrix Poisson ratios
EngConM = C2eng(Cm);
Num     = EngConM(9);

% Eshelby tensors for matrix
EshM = eshtens(aspect, Num);

% Strain concentration tensor 
B = inv4( I4 + EshM*R4*inv4(Cm)*R4*(Cf-Cm) );

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
