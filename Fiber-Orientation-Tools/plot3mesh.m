function plot3mesh(p,elcon,varargin)%PLOTM3MESH Plot a mesh of polygons in 3D.%   PLOT3MESH(P, ELCON) plots the mesh whose nodal coordinates are given%   by P and whose element connectivity is ELCON.%   P(:,I) = the coordinates of node I.%   ELCON(J,K) = the node number I for local node K in element J.%   The number of nodes per element is arbitrary %   (3 is the minimum; 4 is common), but all elements must%   have the same number of nodes.%%   PLOT3MESH(P, ELCON, LINESPEC) plots element boundaries using the %   linespecification of LINESPEC (e.g., 'b-' for solid blue lines).  %%   PLOT3MESH(P, ELCON, SHRINK)reduces the sive of each element by %   the factor SHRINK.  This can be useful in diagnosing element%   connectivity issues.  The default is shrink = 1.%%   The function leaves "hold on" so that subsequent plots%   may be superimposed.%%   See also: PLOTSPHEREMESHwashold = ishold; % Logical; the incoming hold statehold on           % Want hold on during this function % Default values for optional arguments  linespec = 'r-';   shrink = 1;% Handle optional arguments. Note that nargin is the *total* number of args. if (nargin > 2)  for n=1:nargin-2     if ischar(varargin{n})         linespec = varargin{n};     else         shrink = min(varargin{n},1);     end  endend% Find number of elements and number of nodes per element[nelts, nnpe] = size(elcon);for i = 1:nelts   % element centroid coordinates   pc = mean(p(:,elcon(i,:)),2);     pc2 = pc*ones(1,nnpe+1);  % each column of pc2 = pc   % nodal coordinates for shrunken element   pn = pc2 + shrink*(p(:,[elcon(i,:) elcon(i,1)]) - pc2);   % plot the element   plot3(pn(1,:), pn(2,:), pn(3,:), linespec)end% Dress up the plotaxis equalview(3)box ongrid onset(gca, 'FontSize', 16)xlabel('x')ylabel('y')zlabel('z')if ~ washold    % Return the hold state to its previous setting    hold offendreturn