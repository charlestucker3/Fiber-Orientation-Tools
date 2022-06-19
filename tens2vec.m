function Vec = tens2vec(Tens)
%VEC = TENS2VEC(TENS) Convert a symmetric 2nd order tensor TENS, stored in 
%     matrix form, to VEC in column vector (contracted) form.  In 3-D, TENS
%     is 3x3 and VEC is 6x1.  In 2-D, TENS is 2x2 and VEC is 3x1.


% Determine whether Tens is 2-D or 3-D (by its size),
% and contract accordingly
[rows, cols] = size(Tens);
if rows == 3 && cols == 3
    % 3-D tensor
    Vec = [Tens(1,1), Tens(2,2), Tens(3,3), ...
           Tens(2,3), Tens(3,1), Tens(1,2)]';
elseif rows == 2 && cols == 2
    % 2-D tensor
    Vec = [Tens(1,1), Tens(2,2), Tens(1,2)]';
else
    error('Tens must be 2x2 or 3x3')

end