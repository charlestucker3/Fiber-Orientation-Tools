function hs = surfPsi3D(px, py, pz, psi)
%HS = SURFPSI(PX,PY, PZ, PSI) draws a Matlab surf colored surface for the 
%    orientation distribution function PSI on the unit sphere.  Input data
%    is consistent with the finite difference solutions of solvePsi3D, and
%    consists of cell centroid values for a finite volume mesh covering the
%    hemisphere 0 <= theta <= pi, 0 <= phi <= pi.  
%    PX, PY and PZ are the components of the orientation vector at cell
%    centers.  PSI can be a column vector of distribution function values,
%    and will be reshaped to grid form as needed.  
%    The function return HS, a handle to the surf object.

% Get grid dimensions
[ntheta, nphi] = size(px);
if numel(psi) ~= ntheta * nphi
    error('psi must have the same number of elements at px, py, pz')
end

% -- Build the grid data for the full sphere, using periodicity
psi = reshape(psi,ntheta,nphi);  % Half-sphere data, in grid form
psiGrid = [psi, flipud(psi), psi(:,1)]; % Full-sphere data
pxGrid  = [px, -px,  px(:,1)];  
pyGrid  = [py, -py,  py(:,1)];
pzGrid  = [pz,  pz,  pz(:,1)];

% Fill in the "holes" at the poles by adding data rows there
[~, ncol] = size(psiGrid);
psiGrid = [ones(1,ncol)*mean(psiGrid(1,:)); ...
           psiGrid; ...
           ones(1,ncol)*mean(psiGrid(end,:))];
pxGrid  = [zeros(1,ncol); pxGrid; zeros(1,ncol)];
pyGrid  = [zeros(1,ncol); pyGrid; zeros(1,ncol)];
pzGrid  = [ones( 1,ncol); pzGrid; -ones(1,ncol)];

% Make the plot
hs = surf(pxGrid, pyGrid, pzGrid, psiGrid);
axis equal
shading interp

return