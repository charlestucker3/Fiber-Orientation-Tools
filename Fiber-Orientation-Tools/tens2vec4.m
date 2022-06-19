function v4  = tens2vec4(A4)
%V4 = TENS2VEC4(A4) returns a column vector V4 that contains the
%    independent components of a fully symmetric fourth-order tensor A4.
%    In 3-D, A4 is 6x6 and V4 is 15x1.  In 2-D, A4 is 3x3 and V4 is 5x1.
%    While this function is usually applied to orientation tensors, the
%    normalization property of the fourth-order orientation tensor are not
%    used.  However, the input tensor must be fully symmetric.    
%
%    The 3-D tensor components in V4 appear in the following order:
%       Tensor Indices     Matrix Indices
%            1111                11
%            2222                22
%            3333                33
%            2233                23
%            3311                31
%            1122                12
%            1123                14
%            1131                15
%            1112                16
%            2223                24
%            2231                25
%            2212                26
%            3323                34
%            3331                35
%            3312                36
%
%    For 2-D, the tensor components in V4 appear in the following order:
%       Tensor Indices     Matrix Indices
%            1111                11
%            2222                22
%            1122                12
%            1112                13
%            2221                23
%
%  See also: VEC2TENS4

% Determine whether A4 is 2-D or 3-D (by its size), and form v4 accordinly 
[rows, cols] = size(A4);
if rows == 6 && cols == 6
    % 3-D tensor
    v4 = [A4(1,1); A4(2,2); A4(3,3); ...
          A4(2,3); A4(3,1); A4(1,2); ...
          A4(1,4); A4(1,5); A4(1,6); ...
          A4(2,4); A4(2,5); A4(2,6); ...
          A4(3,4); A4(3,5); A4(3,6)];
elseif rows == 3 && cols == 3
    % 2-D tensor
    v4 = [A4(1,1); A4(2,2); A4(1,2); A4(1,3); A4(2,3)];
else
    error('A4 must be 3x3 or 6x6')
end



end