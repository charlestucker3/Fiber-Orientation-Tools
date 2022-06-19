function A4rot = rotate4(A4, Q)
%A4ROT = ROTATE4(A4, Q) gives the fourth-order tensor A4 in a
%        coordinate system rotated by the matrix Q from the original
%        system.  In 3-D, both A4 and A4ROT are stored in contracted form
%        as 6x6 matrices and Q is 3x3.  In 2-D problems A4 and A4ROT are
%        3x3 and Q can be 2x2 or (if the 1-2 plane is the plane of the
%        tensor) 3x3.


% Determine whether A4 is 2-D or 3-D (by its size), and set up the
% appropriate rotation matrix.  
[rows, cols] = size(A4);
if rows == 6 && cols == 6
    % 3-D tensor
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
    
elseif rows == 3 && cols == 3
    % 2-D tensor
    % Form contracted-notation rotation matrix
    Q4a=[Q(1,1)*Q(1,1), Q(1,2)*Q(1,2), Q(1,1)*Q(1,2);
         Q(2,1)*Q(2,1), Q(2,2)*Q(2,2), Q(2,1)*Q(2,2);
         Q(1,1)*Q(2,1), Q(1,2)*Q(2,2), Q(1,1)*Q(2,2)];
    
    Q4b=[0, 0, Q(1,2)*Q(1,1);
         0, 0, Q(2,2)*Q(2,1);
         0, 0, Q(1,2)*Q(2,1)];
    
    Q4  = Q4a+Q4b;  % rotation matrix for   symmetric tensors
else
    error('A4 must be 2x2 or 3x3')
end


% Perform the actual rotation
A4rot = Q4 * A4 * Q4'; % From Nadeau and Ferrari, 1998.

end