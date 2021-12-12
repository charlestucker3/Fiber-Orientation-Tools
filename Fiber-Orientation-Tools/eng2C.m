function C = eng2C(EngCon)
%ENG2C    Convert engineering constants to stiffness tensor.
%    C = ENG2C(ENGCON) converts a set of engineering constants (technical
%    constants) for an orthotropic material to a stiffness tensor C.
%    C is stored in 6x6 matrix form, and the vector of engineering
%    constants is
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

%     Reference: R. M. Jones, Mechanics of Composite Materials,
%     McGraw-Hill, New York (1975), pp40-41.
      
% - Fill in S, the matrix inverse of the stiffness matrix C.  Note that S
%   here is NOT the contraction of the tensor inverse of C.
S = zeros(6,6);

% Tensile and shear moduli are on the diagonal
S(1,1) = 1/EngCon(1);
S(2,2) = 1/EngCon(2);
S(3,3) = 1/EngCon(3);
S(4,4) = 1/EngCon(4);
S(5,5) = 1/EngCon(5);
S(6,6) = 1/EngCon(6);

% Poisson terms are diagonally symmetric
S(1,2) = -1*EngCon(9)/EngCon(1);
S(2,1) =  S(1,2);
S(1,3) = -1*EngCon(8)/EngCon(3);
S(3,1) =  S(1,3);
S(2,3) = -1*EngCon(7)/EngCon(2);
S(3,2) =  S(2,3);
%
% - Invert S to get the C matrix
C = inv(S);
%
return;
