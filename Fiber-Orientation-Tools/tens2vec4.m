function v4  = tens2vec4(A4)
%V4 = TENS2VEC4(A4) returns a 15x1 column vector V4 that contains the
%    independent components of a fully symmetric fourth-order tensor A4.
%    While this function is usually applied to orientation tensors, the
%    normalization property of the fourth-order orientation tensor are not
%    used.  
%
%    The tensor components in V4 appear in the following order:
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
%  See also: VEC2TENS4

v4 = [A4(1,1); A4(2,2); A4(3,3); ...
      A4(2,3); A4(3,1); A4(1,2); ...
      A4(1,4); A4(1,5); A4(1,6); ...
      A4(2,4); A4(2,5); A4(2,6); ...
      A4(3,4); A4(3,5); A4(3,6)];

end