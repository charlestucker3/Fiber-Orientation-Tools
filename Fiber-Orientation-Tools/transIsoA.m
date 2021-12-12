function [A2, A4]= transIsoA(A11, A1111)
%TRANSISOA    Find transversely isotropic orientation tensors
%    [A2, A4] = TRANSISOA(A11, A1111) finds the full second-order
%        orientation tensor A2 and fourth-order orientation tensor A4 for
%        transversely isotropic orientation, using the tensor components
%        A11 and A1111 associated with the symmetry axis.  This relies only
%        on symmetry and normalization.  A2 is returned in 3x3 matrix form
%        and A4 is returned in 6x6 form.

% Second-order tensor
A22 = (1-A11)/2;
A2 = diag([A11, A22, A22]);

% Independent components of the fourth-order tensor
A1122 = (A11 - A1111)/2;
A2233 = (A22 - A1122)/4;
A2222 = 3*A2233;

% Fill in fourth-order tensor using symmetry
A4 = diag([A1111, A2222, A2222, A2233, A1122, A1122]);
A4(1,2) = A1122;
A4(1,3) = A1122;
A4(2,1) = A1122;
A4(2,3) = A2233;
A4(3,1) = A1122;
A4(3,2) = A2233;

return

