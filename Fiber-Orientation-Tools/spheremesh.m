function [p,tricon] = spheremesh(nlevels, varargin)
%[P TRICON] = SPHEREMESH(NLEVELS) generates a triangular mesh on 
%   the unit sphere, by refining an initial octahedral mesh.
%   NLEVELS is the number of refinement levels (0 for no refinement).
%   P (3 x numnod) returns the nodal coordinates and
%   TRICON (numel x 3) returns the node numbers for each element.  
%
%[P TRICON] = SPHEREMESH(NLEVELS, SHAPE) uses the SHAPE as the initial
%   mesh geometry.  Options for SHAPE are:
%       'oct'   Octahedron (the default)
%       'hemi'  Half an octahedron, creating a hemispherical mesh
%       'tet'   Tetrahedron
%       'ico'   Icosahedron (20 triangular faces)
%       'cap'   10 faces, comprising half of the icosahedral mesh
%               This shape has p <-> -p symmetry.  
%
% See also: REFINEMESH, MESHCON

% Programmer's note: The edgecon and triedge matrices for the initial
% (level 0) meshes are hard-coded here, but can be generated for an
% arbitrary mesh using meshcon.m.  

%    The numbers of nodes and elements in the full octahedral mesh are:
%    NLEVELS  Nodes   Elements
%        0        6          8
%        1       18         32
%        2       66        128
%        3      258        512
%        4     1026       2048
%        5     4098       8192
%        6    16386      32768
%        7    65538     131072
%        8   262146     524288
%        9  1048578    2097152

% Set the shape of the initial mesh
if nargin >= 2
    shape = lower(varargin{1});
    shape = shape(1:3);
else
    shape = 'oct'; % The default shape: octahedral
end

% -- Generate the initial mesh (these are all hard-coded)
switch shape
    
    case 'hem'
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

    case 'oct'
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

        
        
% Refine the mesh as many levels as desired
for n = 1:nlevels
    [p, tricon, edgecon, triedge] = refinemesh(p, tricon, edgecon, triedge);
end

        
return;