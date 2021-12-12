function p = sph2p(theta, phi)
%P = SPH2P(THETA, PHI) gives the unit vector P (3x1) corresponding to 
%    spherical-coordinate angles THETA and PHI, which are in radians.
%    If THETA and PHI are vectors of length n, P is 3xn with each column
%    containing a P vector.

np = numel(theta);
p = zeros(3,np);
for i = 1:np
    p(1,i) = cos(phi(i)) * sin(theta(i));
    p(2,i) = sin(phi(i)) * sin(theta(i));
    p(3,i) = cos(theta(i));
end

return