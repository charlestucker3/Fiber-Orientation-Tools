function Aderiv = AdotJeffQuad(Av, L, xi)
%   Aderiv = AdotJeffQuad(Av,L,xi) returns the derivatives Aderiv of the
%   components of the second-order orientation tensor when all fibers
%   follow Jeffery's equation and the quadratic closure is used.  Av
%   contains the orientation tensor components in a column vector, using
%   contracted notation.  L is the velocity gradient tensor (3x3) and xi is
%   the particle shape factor. Aderiv is a 6x1 column vector, also in
%   contracted notation.

A  = vec2tens(Av);           % A in 3x3 matrix form.
A4 = Av*Av';                 % Quadratic closure for the 4th-order 
                             % orientation tensor, (6x6) contracted form
R4 = diag([1 1 1 2 2 2]);    % Used for contracted-notation products
W  = (L-L')/2;               % Vorticity tensor
D  = (L+L')/2;               % Rate-of-deformation tensor

% Main equation
Aderiv = tens2vec(W*A - A*W + xi*(D*A + A*D)) - xi*2*A4*R4*tens2vec(D);

return

