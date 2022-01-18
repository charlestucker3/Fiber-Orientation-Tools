function A4rot = rotate4(A4, Q)
%A4ROT = ROTATE4(A4, Q) gives the fourth-order tensor A4 in a
%        coordinate system rotated by the 3x3 matrix Q from the original
%        system.  Both A4 and A4ROT are stored in contracted form as 6x6
%        matrices.

% Form contracted-notation rotation matrix
Q4a=[Q(1,1)*Q(1,1), Q(1,2)*Q(1,2), Q(1,3)*Q(1,3), ...
     Q(1,2)*Q(1,3), Q(1,3)*Q(1,1), Q(1,1)*Q(1,2);
     Q(2,1)*Q(2,1), Q(2,2)*Q(2,2), Q(2,3)*Q(2,3), ...
     Q(2,2)*Q(2,3), Q(2,3)*Q(2,1), Q(2,1)*Q(2,2);
     Q(3,1)*Q(3,1), Q(3,2)*Q(3,2), Q(3,3)*Q(3,3), ...
     Q(3,2)*Q(3,3), Q(3,3)*Q(3,1), Q(3,1)*Q(3,2);
     Q(2,1)*Q(3,1), Q(2,2)*Q(3,2), Q(2,3)*Q(3,3), ...
     Q(2,2)*Q(3,3), Q(2,3)*Q(3,1), Q(2,1)*Q(3,2);
     Q(3,1)*Q(1,1), Q(3,2)*Q(1,2), Q(3,3)*Q(1,3), ...
     Q(3,2)*Q(1,3), Q(3,3)*Q(1,1), Q(3,1)*Q(1,2);
     Q(1,1)*Q(2,1), Q(1,2)*Q(2,2), Q(1,3)*Q(2,3), ...
     Q(1,2)*Q(2,3), Q(1,3)*Q(2,1), Q(1,1)*Q(2,2)];

Q4b=[0, 0, 0, Q(1,3)*Q(1,2), Q(1,1)*Q(1,3), Q(1,2)*Q(1,1);
     0, 0, 0, Q(2,3)*Q(2,2), Q(2,1)*Q(2,3), Q(2,2)*Q(2,1);
     0, 0, 0, Q(3,3)*Q(3,2), Q(3,1)*Q(3,3), Q(3,2)*Q(3,1);
     0, 0, 0, Q(2,3)*Q(3,2), Q(2,1)*Q(3,3), Q(2,2)*Q(3,1);
     0, 0, 0, Q(3,3)*Q(1,2), Q(3,1)*Q(1,3), Q(3,2)*Q(1,1);
     0, 0, 0, Q(1,3)*Q(2,2), Q(1,1)*Q(2,3), Q(1,2)*Q(2,1)];

Q4  = Q4a+Q4b;  % rotation matrix for   symmetric tensors

% Perform the actual rotation
A4rot = Q4 * A4 * Q4'; % From Nadeau and Ferrari, 1998.

end