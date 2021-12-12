function C = iso2C(E, Nu)
%ISO2C     Convert isotropic elastic properties to a stiffness tensor.
%     C = ISO2C(E, NU) finds the stiffness tensor C, in 6x6 matrix form,
%     for an isotropic material with Young's modulus E and Poisson ratio
%     NU. 

% - Find engineering constants in the form on eng2C.m and C2eng.m
EngCon = zeros(9,1);

% - Extensional moduli
EngCon(1) = E;
EngCon(2) = E;
EngCon(3) = E;

%- Shear moduli
EngCon(4) = E/(2*(1+Nu));
EngCon(5) = EngCon(4);
EngCon(6) = EngCon(4);

% - Poisson Ratios
EngCon(7) = Nu;
EngCon(8) = Nu;
EngCon(9) = Nu;

% - Use general tool to get stiffness tensor
C = eng2C(EngCon);

return
