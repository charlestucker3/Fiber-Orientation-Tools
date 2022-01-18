function [p,tricon, varargout] = spheremesh(nseg, varargin)
%[P TRICON] = SPHEREMESH(NSEG) generates a triangular mesh on 
%   the unit sphere by refining an initial polyhedral mesh of triangles. Each edge of
%   the initial polyhedron is divided into NSEG segments, and a mesh of
%   NSEG^2 triangular elements is generated on each initial triangle. P (3
%   x numnod) returns the nodal coordinates and TRICON (numel x 3) returns
%   the node numbers for each element.  The default initial shape is an
%   octahedron.  
%
%   The default mapping for the octahedron is an equal-area mapping (A.
%   Holhos and D. Rosca, Computers and Math with Applications, 67 (2014)
%   1092-1107).  This spreads the points almost evenly, and all elements
%   have the same area.  It is also possible to use a radial mapping from
%   the initial mesh (p = x/|x|), and this is the only option for other
%   initial shapes (see below).
%
%[P, TRICON, PCEN] = SPHEREMESH(NSEG) also returns the coordinates of the
%   element centroids PCEN (3 x numel).  The centroids are computed on the
%   initial polyhedral faces, then mapped to the sphere using the same
%   mapping as the nodes.
%
%[P TRICON] = SPHEREMESH(NSEG, SHAPE) uses SHAPE as the initial
%   polyhedron.  SHAPE can also control the mapping of points octahedral
%   meshes.  Options for SHAPE are:
%     'oct'     Octahedron, equal-area mapping (the default)
%     'hemi'    Half an octahedron, creating a hemispherical mesh;
%               equal-area mapping
%     'octrad'  Octahedron, radial mapping
%     'hemirad' Half an octahedron, radial mapping
%     'tet'     Tetrahedron, radial mapping
%     'ico'     Icosahedron (20 triangular faces), radial mapping
%     'cap'     10 faces, comprising half of the icosahedral mesh; radial
%               mapping. This mesh covers half of the sphere and has no
%               nodes that are antipodal to other nodes.
%
%   A full octahedral mesh has 8*NSEG^2 elements and 4*NSEG^2 + 2 nodes.
%
% See also: REFINEMESH, MESHCON

% Programmer's note: The edgecon and triedge matrices for the initial
% (level 0) meshes are hard-coded here, but can be generated for an
% arbitrary mesh using meshcon.m.  


% Ensure a reasonable value of nseg
if nseg < 1
    error('NSEG must be >= 1')
end

% Set the shape of the initial mesh
if nargin >= 2
    shape = lower(varargin{1});
else
    shape = 'oct'; % The default shape: octahedral, equal-area mapping
end

% -- Generate the initial mesh (these are all hard-coded)
switch shape
    
    case {'hemi', 'hemirad'}
        % Use only the upper (z >= 0) hemisphere of the octahedron
        p = [1 0 0 -1  0;
             0 1 0  0 -1;
             0 0 1  0  0];
        tricon = [1,2,3;  2,4,3; 4,5,3; 5,1,3];
        % edgecon(i,j) = local node j in edge i
        % Order of nodes for each edge does not matter
        edgecon = [3,1; 3,2; 3,4; 3,5;
                   1,2; 2,4; 4,5; 5,1];
        % triedge(i,j) = edge j in element i
        triedge = [ 5, 2,1;  6, 3,2;  7, 4, 3;  8,1,4];
        % Important: order of edges for each element must match
        % order of nodes in tricon.  E.g., edge 1 connects nodes
        % 1&2; edge 2 connects nodes 2&3; edge 3 connects nodes 3&1.
        
        if strcmp(shape, 'hemi')
            % Increase the size of the octahedron to have the same surface
            % area as the unit sphere
            L = sqrt(2*pi) / 3^(1/4); % Edge length of new octahedron
            p = p * L/sqrt(2);        % Scale up the nodal coordinates
        end

    case {'oct', 'octrad'}
        % Create a full spherical mesh based on an octahedron
        p = [1 0 0 -1  0  0;
             0 1 0  0 -1  0;
             0 0 1  0  0 -1];
        tricon = [1,2,3;  2,4,3; 4,5,3; 5,1,3;
                  1,6,2;  2,6,4; 4,6,5; 5,6,1];
        % edgecon(i,j) = local node j in edge i
        % Order of nodes for each edge does not matter
        edgecon = [3,1; 3,2; 3,4; 3,5;
            1,2; 2,4; 4,5; 5,1;
            1,6; 2,6; 4,6; 6,5];
        % triedge(i,j) = edge j in element i
        triedge = [ 5, 2,1;  6, 3,2;  7, 4, 3;  8,1,4;
                    9,10,5; 10,11,6;  11,12,7; 12,9,8];
        % Important: order of edges for each element must match
        % order of nodes in tricon.  E.g., edge 1 connects nodes
        % 1&2; edge 2 connects nodes 2&3; edge 3 connects nodes 3&1.
        
        if strcmp(shape, 'oct')
            % Increase the size of the octahedron to have the same surface
            % area as the unit sphere
            L = sqrt(2*pi) / 3^(1/4); % Edge length of new octahedron
            p = p * L/sqrt(2);        % Scale up the nodal coordinates
        end
        
    case 'ico'
        % Create an icosahedral mesh
        % -- Nodal coords in (theta, phi), as row vectors
        ptheta = [0, (pi/2 - atan(1/2))*ones(1,5), ...
                     (pi/2 + atan(1/2))*ones(1,5), pi];
        pphi   = [0, (0:2:8)*pi/5, (1:2:9)*pi/5, 0];
        
        % -- Convert to Cartesian coords
        p      = zeros(3,12);
        p(1,:) = cos(pphi) .* sin(ptheta);
        p(2,:) = sin(pphi) .* sin(ptheta);
        p(3,:) = cos(ptheta);
        
        % -- Create element connectivity
        tricon = [ 1,  2,  3;
                   1,  3,  4;
                   1,  4,  5;
                   1,  5,  6;
                   1,  6,  2;
                   2,  7,  3;
                   3,  8,  4;
                   4,  9,  5;
                   5, 10,  6;
                   6, 11,  2;
                   7,  8,  3;
                   8,  9,  4;
                   9, 10,  5;
                  10, 11,  6;
                  11,  7,  2;
                  12,  7,  8;
                  12,  8,  9;
                  12,  9, 10;
                  12, 10, 11;
                  12, 11,  7];
        
        % -- Edge connectivity (node numbers for each edge)
        edgecon = [1,  2;
                   1,  3;
                   1,  4;
                   1,  5;
                   1,  6;
                   2,  3;
                   3,  4;
                   4,  5;
                   5,  6;
                   6,  2;
                   2,  7;
                   3,  8;
                   4,  9;
                   5, 10;
                   6, 11;
                   2, 11;
                   3,  7;
                   4,  8;
                   5,  9;
                   6, 10;
                   7,  8;
                   8,  9;
                   9, 10;
                  10, 11;
                  11,  7;
                   7, 12;
                   8, 12;
                   9, 12;
                  10, 12;
                  11, 12];
        
      
        % -- Element edges
        triedge = [ ...
            1     6     2;
            2     7     3;
            3     8     4;
            4     9     5;
            5    10     1;
            11    17     6;
            12    18     7;
            13    19     8;
            14    20     9;
            15    16    10;
            21    12    17;
            22    13    18;
            23    14    19;
            24    15    20;
            25    11    16;
            26    21    27;
            27    22    28;
            28    23    29;
            29    24    30;
            30    25    26];
        
              
    case 'cap'
        % Create the 10-face "cap" of an icosahedral mesh
        % -- Nodal coords in (theta, phi), as row vectors
        %    (Uses all the ico nodes except #12).
        ptheta = [0, (pi/2 - atan(1/2))*ones(1,5), ...
                     (pi/2 + atan(1/2))*ones(1,5)];
        pphi   = [0, (0:2:8)*pi/5, (1:2:9)*pi/5];
        
        % -- Convert to Cartesian coords
        p      = zeros(3,11);
        p(1,:) = cos(pphi) .* sin(ptheta);
        p(2,:) = sin(pphi) .* sin(ptheta);
        p(3,:) = cos(ptheta);
        
        % -- Create element connectivity
        tricon = [ 1,  2,  3;
                   1,  3,  4;
                   1,  4,  5;
                   1,  5,  6;
                   1,  6,  2;
                   2,  7,  3;
                   3,  8,  4;
                   4,  9,  5;
                   5, 10,  6;
                   6, 11,  2];
        
        % -- Edge connectivity (node numbers for each edge)
        edgecon = [1,  2;
                   1,  3;
                   1,  4;
                   1,  5;
                   1,  6;
                   2,  3;
                   3,  4;
                   4,  5;
                   5,  6;
                   6,  2;
                   2,  7;
                   3,  8;
                   4,  9;
                   5, 10;
                   6, 11;
                   2, 11;
                   3,  7;
                   4,  8;
                   5,  9;
                   6, 10];
        
      
        % -- Element edges
        triedge = [ ...
            1     6     2;
            2     7     3;
            3     8     4;
            4     9     5;
            5    10     1;
            11    17     6;
            12    18     7;
            13    19     8;
            14    20     9;
            15    16    10];

    case 'tet'
        p = [  0,           0,        1;  % From Wikipedia, "Tetrahedron"
            sqrt(8/9),      0,     -1/3;
           -sqrt(2/9),  sqrt(2/3), -1/3;
           -sqrt(2/9), -sqrt(2/3), -1/3]';
        tricon = [1, 2, 3;
                  1, 3, 4;
                  1, 4, 2;
                  2, 3, 4];
        edgecon = [1 2; 1, 3; 1, 4; 2, 3; 3, 4; 4, 2];
        triedge = [1, 4, 2;
                   2, 5, 3;
                   3, 6, 1;
                   4, 5, 6];
                  
    otherwise
        error('Illegal value of SHAPE')
end
% ---------- simple mesh for debugging ------
%tricon = [1,2,3; 2,4,3];
%edgecon = [1,2; 2,3; 3,1; 2,4; 4,3];
%triedge = [1,2,3; 4,5,2]; % good triedge data
%triedge = [1,2,3; 5,2,4]; % edge order ~= node order; refinement fails!
% -------------------------------------------

        
        
% Refine the mesh.  If nseg == 1, refinemesh returns the original mesh
% along with the element centroids.
[p, tricon, pcen] = refinemesh(p, tricon, edgecon, triedge, nseg);

% Map the polyhedral mesh to the sphere
switch shape
    case {'oct', 'hemi'}
        % Equal-area mapping.  This version of the Holhos and Rosca map is
        % given by Hardin, Michaels and Saff (2016),
        % https://arxiv.org/abs/1607.04590 and avoids some divide-by-zero
        % issues with the equations from the original paper.
        
        % Nodes
        p(3,:) = 2*p(3,:) .* (sqrt(2)*L - abs(p(3,:))) / L^2;
        gamma = (pi/2) * abs(p(2,:)) ./ (abs(p(1,:)) + abs(p(2,:)));
        gamma(p(1,:)==0 & p(2,:)==0) = 0;  % Fix divide-by-zero @ N and S poles
        p(1,:) = sign(p(1,:)) .* sqrt(1-p(3,:).^2) .* cos(gamma);
        p(2,:) = sign(p(2,:)) .* sqrt(1-p(3,:).^2) .* sin(gamma);
        
        % Centroids
        pcen(3,:) = 2*pcen(3,:) .* (sqrt(2)*L - abs(pcen(3,:))) / L^2;
        gamma = (pi/2) * abs(pcen(2,:)) ./ (abs(pcen(1,:)) + abs(pcen(2,:)));
        gamma(pcen(1,:)==0 & pcen(2,:)==0) = 0;  % Fix divide-by-zero
        pcen(1,:) = sign(pcen(1,:)) .* sqrt(1-pcen(3,:).^2) .* cos(gamma);
        pcen(2,:) = sign(pcen(2,:)) .* sqrt(1-pcen(3,:).^2) .* sin(gamma);
        
    otherwise
        % Radial mapping
        p    = p    ./ vecnorm(p);
        pcen = pcen ./ vecnorm(pcen);
end

% Return the element centroids, if requested
if nargout >= 3
    varargout{1} = pcen;
end

return;