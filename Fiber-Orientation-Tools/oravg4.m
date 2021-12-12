function Cavg = oravg4(C, A4, varargin)
%ORAVG    Orientation average of a fourth-order tensor
%     CAVG = ORAVG4(C, A4) returns CAVG, the orientation average of the
%     transversly isotropic fourth-order tensor C, using the fourth-order
%     orientation tensor A4.  C, A4 and CAVG are stored in 6x6 matrix
%     (contracted) form.
%
%     C may be a NON-SYMMETRIC matrix, e.g. a strain concentration tensor.
%     C and A4 must be referenced to the same axes, and C must be
%     transversely isotropic with 1 as the symmetry axis.
%
%     A4 should have the normalization property of the exact 4th-order 
%     tensor, A_ijkk = A_ij.  This is not true for a few closure
%     approximations, in which case use
%     CAVG = ORAVG(C, A4, A2) where A2 is the second-order orientation
%     tensor in 3x3 matrix form.  
%

%     Reference for symmetric C tensors: S.G. Advani and C.L. Tucker, 
%     J. Rheology, vol. 31, pp751-784 (1987).

%  Kronecker Delta function
delta = eye(3);

%  Tables to expand and contract indices (i.e., translate between
%  Cartesian tensor and contracted notation)
iexpand = [1, 2, 3, 2, 3, 1];
jexpand = [1, 2, 3, 3, 1, 2];
icontract  = [1, 6, 5; 6, 2, 4; 5, 4, 3];

%  Find the 6 invariants of the input C tensor
B1  = C(1,1) + C(2,2) - C(1,2) - C(2,1) - 4*C(6,6);
B2U = C(1,2) - C(2,3);
B2L = C(2,1) - C(2,3);
B3  = C(6,6) + (C(2,3)-C(2,2))/2;
B4  = C(2,3);
B5  = (C(2,2)-C(2,3))/2;

if nargin >= 3
    % The second-order tensor was provided on input
    A2 = varargin{1};
else
    %  Find the components of the second-order tensor A2 in tensor notation
    %  using the normalization property of A4.
    A2 = zeros(3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                A2(i,j) = A2(i,j) + A4(icontract(i,j),icontract(k,k));
            end
        end
    end
end

%  Compute the orientation average component by component.
%  Note: ijkl are the expanded indices for mn.
Cavg = zeros(6);
for m = 1:6
    for n = 1:6
        i = iexpand(m);
        j = jexpand(m);
        k = iexpand(n);
        l = jexpand(n);
        Cavg(m,n) = B1  * A4(m,n) ...
            + B2U * (A2(i,j)*delta(k,l)) ...
            + B2L * (A2(k,l)*delta(i,j)) ...
            + B3  * (A2(i,k)*delta(j,l) + A2(i,l)*delta(j,k)  ...
            + A2(j,l)*delta(i,k) + A2(j,k)*delta(i,l)) ...
            + B4  * delta(i,j)*delta(k,l) ...
            + B5  *(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
    end
end

return

