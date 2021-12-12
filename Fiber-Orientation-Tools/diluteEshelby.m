function C = diluteEshelby(Cf, Cm, vf, aspect)
%DILUTEESHELBY     Dilute Eshelby model for stiffness, aligned fibers.
%      C = DILUTE(CF, CM, VF, ASPECT) uses the dilute form of the Eshelby
%      model to compute the stiffness tensor C of a composite with
%      unidirectional reinforcement. Cf and Cm are the stiffness tensors
%      for the fiber and matrix in 6x6 form, and C is returned in the same
%      form.  VF is the fiber volume fraction and ASPECT is the fiber
%      aspect ratio (length/diameter) or the ellipsoidal aspect ratio.
%
%      Cm must be isotropic and Cf must be either isotropic or transversely
%      anisotropic about the 1 axis.  C will be transversely isotropic
%      about the 1 axis. 
 

I4 = diag([1,1,1,0.5,0.5,0.5]);  % 4th-order unit tensor
R4 = diag([1,1,1,2,  2,  2]);    % Used in contracted-notation products


% Matrix Poisson ratio
EngCon = C2eng(Cm);
Num    = EngCon(9);

% Eshelby tensor
Esh = eshtens(aspect, Num);

% Strain concentration tensor for Eshelby dilute model
H = inv4( I4 + Esh*R4*inv4(Cm)*R4*(Cf-Cm) );

% Composite stiffness
C = Cm + vf*(Cf - Cm) * R4 * H;

return
