function EngCon = C2eng(C)
%C2ENG   Convert stiffness tensor to engineering constants.
%    ENGCON = C2ENG(C) finds a 9x1 vector of engineering constants
%    (technical constants) ENGCON corresponding to the fourth-order
%    elastic stiffness tensor C.  C is stored as a 6x6 matrix.
%
%    The engineering constants are:
%          EngCon(1) = E11
%          EngCon(2) = E22
%          EngCon(3) = E33
%          EngCon(4) = G23
%          EngCon(5) = G31
%          EngCon(6) = G12
%          EngCon(7) = Nu23
%          EngCon(8) = Nu31
%          EngCon(9) = Nu12
%
%    If C is not orthotropic and/or the coordinate axes are not the
%    principal material axes then there can be other non-zero elastic
%    properties not represented in ENGCON.

% - Find C4-inverse numerically.  Note that this is NOT the contraction of
%   the tensor inverse of C.
S = inv(C);

% - Get the engineering constants from the components of S.
EngCon = zeros(9,1);
%
%   Extensional moduli
EngCon(1) = 1/S(1,1);
EngCon(2) = 1/S(2,2);
EngCon(3) = 1/S(3,3);

%   Shear moduli
EngCon(4) = 1/S(4,4);
EngCon(5) = 1/S(5,5);
EngCon(6) = 1/S(6,6);

%  Poisson Ratios
EngCon(7) = -S(2,3)*EngCon(2);
EngCon(8) = -S(3,1)*EngCon(3);
EngCon(9) = -S(1,2)*EngCon(1);

return
