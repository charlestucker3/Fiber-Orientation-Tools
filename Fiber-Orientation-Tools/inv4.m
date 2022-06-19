function S = inv4(C)
%INV4    Inverse of a fourth-order tensor in matrix form.
%    S = INV4(C) returns the tensor inverse S of a fourth-order tensor C,
%    with both tensors in contracted form.  In 3-D, both S and C are 6x6
%    matrices; in 2-D they are 3x3.  Note that INV4(C) is NOT simply the
%    matrix inverse of C.

% Determine whether C is 2-D or 3-D (by its size), and set up the
% appropriate symmetric identity tensor I4.  
[rows, cols] = size(C);
if rows == 6 && cols == 6
    % 3-D tensor
    I4 = diag([1,1,1,0.5,0.5,0.5]);
elseif rows == 3 && cols == 3
    % 2-D tensor
    I4 = diag([1,1,0.5]);
else
    error('C must be 2x2 or 3x3')
end

% Find the tensor inverse.
S = I4 * inv(C) * I4;

return
