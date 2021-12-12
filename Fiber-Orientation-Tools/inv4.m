function S = inv4(C)
%INV4    Inverse of a fourth-order tensor, 6x6 matrix form.
%    S = INV4(C) returns the tensor inverse S of a fourth-order tensor
%    C.  Both S and C are stored as 6x6 matrices in contracted form.
%    Note that INV4(C) is NOT simply the matrix inverse of C.

% I4 accounts for the fact that 6x6 matrix multiplication does not count
% the terms the same number of times as indicial notation for fourth-order
% tensors.
I4 = diag([1,1,1,0.5,0.5,0.5]);
S = I4 * inv(C) * I4;

return
