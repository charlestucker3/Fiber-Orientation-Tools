function Tens = vec2tens(Vec)
% TENS = VEC2TENS(VEC) converts a symmetric tensor VEC from column-vector
%       (contracted) form to a tensor form TENS.  In 3-D, VEC is 6x1 and
%       TENS is 3x3.  In 2-D, VEC is 3x1 and TENS is 2x2.

% Determine whether Vec is 2-D or 3-D (by its size),
% and expand accordingly
[rows, cols] = size(Vec);
if (rows == 6 && cols == 1) || (rows == 1 && cols == 6) 
    % 3-D tensor
    Tens = [Vec(1), Vec(6), Vec(5);
           Vec(6), Vec(2), Vec(4);
           Vec(5), Vec(4), Vec(3)];
elseif (rows == 3 && cols == 1) || (rows == 1 && cols == 3) 
    % 2-D tensor
    Tens = [Vec(1), Vec(3);
            Vec(3), Vec(2)];
else
    error('Vec must be 6x1 or 3x1')

end

end