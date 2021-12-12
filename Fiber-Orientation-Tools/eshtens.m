function E = eshtens(alpha, Nu)
%ESHTENS    Eshelby tensor.
%     E = ESHTENS(ALPHA, NU) calculatse Eshelby's tensor E for a spheroid
%     having aspect ratio ALPHA and in an isotropic matrix with Poisson
%     ratio NU.  E is returned in contracted notation as a 6x6 matrix, with
%     the 1 axis is the axis of symmetry.
%

%     Reference: 
%     M.Taya and R. J. Arsenault, "Metal Matrix Composites: 
%     Thermomechanical Behavior," Pergamon Press, pp250-254 (1989).

%  Set elements of E not otherwise calculated to zero
E = zeros(6);
%
%  Sphere, prolate, and oblate cases must be handled separately
%  The parameter "eps" sets the level for activating the sphere formulas
eps = 0.0001;

if (abs(alpha - 1) < eps)
    %        -- Spherical particle
    E(1,1) = (7 - 5*Nu) / (15*(1-Nu));
    E(2,2) = E(1,1);
    E(3,3) = E(1,1);
    E(1,2) = (5*Nu-1) / (15*(1-Nu));
    E(2,1) = E(1,2);
    E(2,3) = E(1,2);
    E(3,2) = E(1,2);
    E(3,1) = E(1,2);
    E(1,3) = E(1,2);
    E(6,6) = (4 - 5*Nu) / (15*(1-Nu));
    E(4,4) = E(6,6);
    E(5,5) = E(6,6);
    
elseif (alpha > 1)
    % -- Prolate (fiber-like) ellipsoid
    G = alpha * ( alpha*sqrt(alpha^2-1) - acosh(alpha) ) / (alpha^2-1)^(1.5);
    
    E(1,1) = (   1 - 2*Nu + (3*alpha^2-1)/(alpha^2-1) ...
        - (1 - 2*Nu + (3*alpha^2  )/(alpha^2-1)) * G ) / (2*(1-Nu));
    E(2,2) = (3*alpha^2) / (8*(1-Nu)*(alpha^2-1)) ...
        + ( 1 - 2*Nu - 9/(4*(alpha^2-1)) ) * G / (4*(1-Nu));
    E(3,3) = E(2,2);
    E(2,3) = (alpha^2/(2*(alpha^2-1)) - G*(1 - 2*Nu + 3/(4*(alpha^2-1)))) ...
        / (4*(1-Nu));
    E(3,2) = E(2,3);
    E(2,1) = - alpha^2/(2*(1-Nu)*(alpha^2-1)) ...
        + ( 3*alpha^2/(alpha^2-1) - (1-2*Nu) ) * G / (4*(1-Nu));
    E(3,1) = E(2,1);
    E(1,2) = - (1 - 2*Nu + 1/   (alpha^2-1) )     / (2*(1-Nu)) ...
        + (1 - 2*Nu + 3/(2*(alpha^2-1))) * G / (2*(1-Nu));
    E(1,3) = E(1,2);
    E(4,4) =  alpha^2 / (8*(1-Nu)*(alpha^2-1)) ...
        + (1 - 2*Nu - 3/(4*(alpha^2-1))) * G / (4*(1-Nu));
    E(5,5) =  ( 1 - 2*Nu -   (alpha^2+1)/(alpha^2-1) )     / (4*(1-Nu)) ...
        - ( 1 - 2*Nu - 3*(alpha^2+1)/(alpha^2-1) ) * G / (8*(1-Nu));
    E(6,6) = E(5,5);
    
else
    % -- Oblate (disk-like) ellipsoid
    G = alpha * ( acos(alpha) - alpha*sqrt(1-alpha^2) ) / (1-alpha^2)^(1.5);
    
    E(1,1) =   ( 4 - 2*Nu - 2/(1-alpha^2))     / (2*(1-Nu)) ...
        + (-4 + 2*Nu + 3/(1-alpha^2)) * G / (2*(1-Nu));
    E(2,2) = - (3*alpha^2) / (8*(1-Nu)*(1-alpha^2)) ...
        + ( 1 - 2*Nu + 9/(4*(1-alpha^2)) ) * G / (4*(1-Nu));
    E(3,3) = E(2,2);
    E(2,3) =  (1 - 1/(1-alpha^2)) / (8*(1-Nu)) ...
        + (-4*(1-2*Nu) + 3/(1-alpha^2) ) * G / (16*(1-Nu));
    E(3,2) = E(2,3);
    E(2,1) =  alpha^2/(2*(1-Nu)*(1-alpha^2)) ...
        - ( (1-2*Nu) +  3*alpha^2/(1-alpha^2) ) * G / (4*(1-Nu));
    E(3,1) = E(2,1);
    E(1,2) =   ( -(1-2*Nu) + 1/(1-alpha^2))     / (2*(1-Nu)) ...
        + (2*(1-2*Nu) - 3/(1-alpha^2)) * G / (4*(1-Nu));
    E(1,3) = E(1,2);
    E(4,4) =  - alpha^2 / (8*(1-Nu)*(1-alpha^2)) ...
        + ( 3/(1-alpha^2) + 4*(1-2*Nu) ) * G / (16*(1-Nu));
    E(5,5) =  ( 1 - 2*Nu +   (1+alpha^2)/(1-alpha^2) )     / (4*(1-Nu)) ...
        - ( 1 - 2*Nu + 3*(1+alpha^2)/(1-alpha^2) ) * G / (8*(1-Nu));
    E(6,6) = E(5,5);
end

return
