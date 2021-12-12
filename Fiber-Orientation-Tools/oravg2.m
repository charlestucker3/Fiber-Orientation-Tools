function Tavg = oravg2(T2, A2)
%ORAVG2    Orientation average of a second-order tensor
%     TAVG = ORAVG2(T2, A2) returns TAVG, the orientation average of the
%     transversly isotropic second-order tensor T, using the second-order
%     orientation tensor A2.  T2 and TAVG are stored in 6x1 column vector
%     form (contracted notation) while A2 is stored in 3x3 matrix form. T2
%     must have the 1 axis as its axis of symmetry.  This implies that T2
%     must be diagonal with the second and third diagonal elements equal.
%

% Check that T2 is transversely isotropic with 1 as its symmetry axis
if sum(abs(T2(4:6))) > 1e-12
    error('T2 must be diagonal')
end
if abs(T2(2) - T2(3))/max(abs(T2(2:3))) > 1e-12
    error('T2 must be transversely isotropic with 1 as its symmetry axis')
end

% Invariants of T2
a1 = T2(1) - T2(2);
a2 = T2(2);

% Tavg = a1 * A2 + a2 * I2, where I2 is the identity tensor
Tavg = a1*tens2vec(A2) + a2 * [1,1,1,0,0,0]';

return

