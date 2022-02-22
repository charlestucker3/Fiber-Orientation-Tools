function varargout = drawPsi(F)
%DRAWPSI(F) draws a sphere colored by the orientation distribution
%    function psi(p) corresponding to the deformation gradient tensor F
%    (3x3), using the Jeffery distribution function.  Use A2F to find F
%    from an orienation tensor A.
%
%HS = DRAWPSI(F) returns a handle to the surf object of the sphere.

% -- Build a toppologically rectangular grid on the unit sphere
nsph = 80;  % Density of sphere grid
[x, y, z] = sphere(nsph);
        

% -- Find the orientation distribution function psi at each point
if cond(F) < 1e6
    Binv = inv(F*F');  % Binv is Montgomery-Smith's B tensor
else  % Avoid warning when F and B are ill-conditioned
    [R, Evals] = eig(F*F');
    Binv = R * diag(1./diag(Evals)) * R';
end
psi = zeros(size(x));
for i = 1:nsph+1
    for j = 1:nsph+1
        p = [x(i,j); y(i,j); z(i,j)];  % p vector for this point
        psi(i,j) = 1 / (4*pi * (p' * Binv * p)^(3/2));
    end
end

% -- Draw the sphere, colored by the distribution function
hs = surf(x, y, z, psi);
shading interp
axis equal

if nargout >= 1  % Return the handle to the surf object
    varargout{1} = hs;
end

return