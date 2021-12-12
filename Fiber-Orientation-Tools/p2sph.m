function [theta, phi] = p2sph(p)
%[THETA, PHI] = P2SPH(P) finds the spherical-coordinate angles THETA and
%     PHI corresponding to the vector P, which is 3x1.  THETA and PHI are
%     returned in radians.
%     If P is 3xn, THETA and PHI are 1xn.

theta = acos(p(3,:));
phi   = atan2(p(2,:), p(1,:));

return