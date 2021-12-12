function [A2, varargout] = F2A(F, varargin)
%A2 = F2A(F) computes the second order orientation
%     tensor A2 (3x3) by applying the deformation gradient tensor F to a
%     uniform initial distribution, using the exact scheme of
%     Montgomery-Smith  (S. Montgomery-Smith, S.,W. He, D. A. Jack and D.
%     S. Smith, "Exact tensor closures for the three-dimensional Jeffery's
%     equation."J. Fluid Mechanics, 680, 312-335, (2011).
%
%[A2, A4] = F2A(F) also computes the fourth-order tensor A4, in 6x6 form.
%     
%     The package elfun18 must be installed, to provide Carlson elliptic
%     integrals. 
%     Milan Batista (2020). elfun18,
%     https://www.mathworks.com/matlabcentral/fileexchange/65915-elfun18,
%     MATLAB Central File Exchange. Retrieved July 17, 2020.

%[A2, A4] = F2A(F, BTOL) is an undocumented option that also sets the
%     value of BTOL, which activates the special A4 approximations when two
%     or more eigenvalues are nearly equal (and avoids divide-by- zero
%     errors).  A2 is unaffected by BTOL. The default value, BTOL = 1e-4,
%     was tuned to give very small errors in A2 (~1e-9) at the changeover
%     to the special formulae (where the error is largest) without exciting
%     rounding errors in the main formula.

% -- Protect against NaN or Inf values of F
Fnan = isnan(F);
if sum(sum(Fnan)) ~= 0
    fprintf('NaN F in F2A.m:\n')
    fprintf('   %f   %f   %f\n', F');
    error('Cannot process NaN F')
end
Finf = isinf(F);
if sum(sum(Finf)) ~= 0
    fprintf('Inf F in F2A.m:\n')
    fprintf('   %f   %f   %f\n', F');
    error('Cannot process Inf F')
end

% -- Find the Finger tensor B.  This is the inverse of Mongomery-Smith's
%    B, but we'll invert it once we have the eigenvalues
B = F*F';

% --Find principal values and axes of B (as in ortho.m)
[Evecs, Bprin] = eig(B);    % eigenvectors/values of B
% b(1), b(2), b(3) are the eigenvalues of my inv(B), and
% match b1, b2, b3 in Montgomery-Smith.  Sort them in 
% increasing order, which will produce the eigenvalues of A2 in 
% decreasing order.  
[b, sortOrder] = sort(1./diag(Bprin)); 
% Normalize the b(i) values so that F, B, etc. have determinant = 1. The
% input value of det(F) can be a bit off due to numerical errors, and not
% having this = 1 on input can alter the output A2 substantially,
% even though the physics of orientation says that det(F) has no effect.
b = b/prod(b)^(1/3);
% Sorted eigenvectors form the rotation matrix
R = [Evecs(:,sortOrder(1)), Evecs(:,sortOrder(2)), Evecs(:,sortOrder(3))];
% Now B = R * diag(1./b) * R'
               

% -- The eigenvalues of A are related to b by the R_D form of Carlson
%    elliptic integrals.  
lam1 = CarlsonRD(b(2), b(3), b(1))/3;
lam2 = CarlsonRD(b(3), b(1), b(2))/3;
lam3 = CarlsonRD(b(1), b(2), b(3))/3;

% -- Return A2 to the original coordinate system
A2 = R * diag([lam1, lam2, lam3]) * R';

% ----- That's all for A2.  The rest of the function deals with A4 -----

if nargout >= 2   % Then also calculate the fourth-order tensor A4
    A4 = zeros(6,6);  % We'll build this initially in principal axes.
    
    % -- First, fill in some main components
    %    Some special forms are needed if any two of the b(i) values 
    %    are very close to one another. 
    if nargin >= 2 % Apply the user value of btol
        btol = varargin{1};
    else % Use the default value
        btol = 1e-4;  % Tolerance for eigenvalues of b being nearly equal.
        % Approximate formulas have errors ~ (btol^2)/10.
    end

    if b(3) - b(1) < btol  % Then all 3 eigenvalues are very close
        c = b - 1;  % To match c1, c2, c3 in Montgomery-Smith
%         fprintf('F2A: 3 eigenvalues nearly equal\n')
        A4(1,1) = 1/5  - 3*c(1)/14 - 3*c(2)/70 - 3*c(3)/70;
        A4(2,2) = 1/5  - 3*c(1)/70 - 3*c(2)/14 - 3*c(3)/70;
        A4(3,3) = 1/5  - 3*c(1)/70 - 3*c(2)/70 - 3*c(3)/14;
        A4(1,2) = 1/15 - 3*c(1)/70 - 3*c(2)/70 - 1*c(3)/70;
        A4(1,3) = 1/15 - 3*c(1)/70 - 1*c(2)/70 - 3*c(3)/70;
        A4(2,3) = 1/15 - 1*c(1)/70 - 3*c(2)/70 - 3*c(3)/70;
        
    elseif b(3) - b(2) < btol  % Then b(2) and b(3) are close
%         fprintf('F2A: b(2) nearly equal b(3)\n')
        A4(1,2) = (b(1)*lam1 - b(2)*lam2) / (2*(b(1)-b(2)));
        A4(1,3) = (b(1)*lam1 - b(3)*lam3) / (2*(b(1)-b(3)));
        A4(2,3) = (lam2 + lam3 - A4(1,2) - A4(1,3))/8;
        A4(1,1) = lam1 - A4(1,2) - A4(1,3);
        A4(2,2) = lam2 - A4(1,2) - A4(2,3);
        A4(3,3) = lam3 - A4(1,3) - A4(2,3);
        
    elseif b(2) - b(1) < btol  % Then b(1) and b(2) are close
%         fprintf('F2A: b(1) nearly equal b(2)\n')
        A4(1,3) = (b(1)*lam1 - b(3)*lam3) / (2*(b(1)-b(3)));
        A4(2,3) = (b(2)*lam2 - b(3)*lam3) / (2*(b(2)-b(3)));
        A4(1,2) = (lam1 + lam2 - A4(1,3) - A4(2,3))/8;
        A4(1,1) = lam1 - A4(1,2) - A4(1,3);
        A4(2,2) = lam2 - A4(1,2) - A4(2,3);
        A4(3,3) = lam3 - A4(1,3) - A4(2,3);
        
    else % The main case, when all the b(i) are distinct
%         fprintf('F2A: all b(i) values are distinct\n')
        A4(1,2) = (b(1)*lam1 - b(2)*lam2) / (2*(b(1)-b(2)));
        A4(2,3) = (b(2)*lam2 - b(3)*lam3) / (2*(b(2)-b(3)));
        A4(1,3) = (b(1)*lam1 - b(3)*lam3) / (2*(b(1)-b(3)));
        A4(1,1) = lam1 - A4(1,2) - A4(1,3);
        A4(2,2) = lam2 - A4(1,2) - A4(2,3);
        A4(3,3) = lam3 - A4(1,3) - A4(2,3);
        
    end
    
    % -- Fill in the rest of A4 by symmetry
    A4(2,1) = A4(1,2);
    A4(3,1) = A4(1,3);
    A4(3,2) = A4(2,3);
    A4(4,4) = A4(2,3);
    A4(5,5) = A4(3,1);
    A4(6,6) = A4(1,2);
    
    % -- Now rotate A4 to the initial axis system (as in ortho.m)
    Id4 = diag( [1 1 1 0.5 0.5 0.5], 0);  % 4th-order identity tensor
    R4  = diag( [1 1 1 2   2   2  ], 0);  % inverse of Id4
    %---- 6x6 rotation matrices and associated quantities -----
    Qa=[R(1,1)*R(1,1), R(1,2)*R(1,2), R(1,3)*R(1,3), ...
        R(1,2)*R(1,3), R(1,3)*R(1,1), R(1,1)*R(1,2);
        R(2,1)*R(2,1), R(2,2)*R(2,2), R(2,3)*R(2,3), ...
        R(2,2)*R(2,3), R(2,3)*R(2,1), R(2,1)*R(2,2);
        R(3,1)*R(3,1), R(3,2)*R(3,2), R(3,3)*R(3,3), ...
        R(3,2)*R(3,3), R(3,3)*R(3,1), R(3,1)*R(3,2);
        R(2,1)*R(3,1), R(2,2)*R(3,2), R(2,3)*R(3,3), ...
        R(2,2)*R(3,3), R(2,3)*R(3,1), R(2,1)*R(3,2);
        R(3,1)*R(1,1), R(3,2)*R(1,2), R(3,3)*R(1,3), ...
        R(3,2)*R(1,3), R(3,3)*R(1,1), R(3,1)*R(1,2);
        R(1,1)*R(2,1), R(1,2)*R(2,2), R(1,3)*R(2,3), ...
        R(1,2)*R(2,3), R(1,3)*R(2,1), R(1,1)*R(2,2)];
    
    Qb=[0, 0, 0, R(1,3)*R(1,2), R(1,1)*R(1,3), R(1,2)*R(1,1);
        0, 0, 0, R(2,3)*R(2,2), R(2,1)*R(2,3), R(2,2)*R(2,1);
        0, 0, 0, R(3,3)*R(3,2), R(3,1)*R(3,3), R(3,2)*R(3,1);
        0, 0, 0, R(2,3)*R(3,2), R(2,1)*R(3,3), R(2,2)*R(3,1);
        0, 0, 0, R(3,3)*R(1,2), R(3,1)*R(1,3), R(3,2)*R(1,1);
        0, 0, 0, R(1,3)*R(2,2), R(1,1)*R(2,3), R(1,2)*R(2,1)];
    
    Q  = Qa+Qb;  % rotation matrix for   symmetric tensors
    
    A4 = Q *A4*R4*inv(Q)*Id4;  % Fourth-order tensor in initial coords.
    varargout{1} = A4;
end

return