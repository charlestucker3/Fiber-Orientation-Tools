function [pNew, triconNew, edgeconNew, triedgeNew] = ...
             refinemesh(pOld, triconOld, edgeconOld, triedgeOld, varargin)
%[PNEW, TRICONNEW, EDGECONNEW, TRIEDGENEW] = REFINEMESH(POLD, TRICONOLD,
%   EDGECONOLD, TRIEDGEOLD) refines a triangular mesh on 3-D unit sphere
%   by adding a node at the midpoint of each edge and creating four new
%   triangles from each old one. 
%   P (3xn) contains the nodal coordinates, TRICON (nelt x 3) lists the
%   nodes for each element, EDGECON (nedge x 2) lists the nodes for each
%   edge and TRIEDGE (nelt x 3) lists the edges in each element.  The 
%   order of edges in TRIEDGE must match the order of nodes in TRICON:
%   Edge 1 connects local nodes 1 and 2, edge 2 connects local nodes 
%   2 and 3, and edge 3 connects local nodes 3 and 1.
%
%   EDGECON and TRIEDGE can be generated from TRICON using the
%   MESHCON function.
%
%   If a fifth argment is present it gives the SHAPE of the mesh. The
%   default is 'sphere', which causes all the PNEW vectors to be normalized
%   to unit length.  A value of 'flat' will keep each new element in the
%   plane of the old element so, for example, a planar mesh remains planar.

% -- Set a logical flag for the shape of the mesh
if nargin >= 5
    if strcmpi(varargin{1}, 'flat')
        sphere = false;
    elseif strcmpi(varargin{1}, 'sph')
        sphere = true;
    else
        error('Illegal value of SHAPE')
    end
else
    sphere = true;  % The default
end

% Work out old and new sizes, and create new storage space
[~,nnold]    = size(pOld);       % nnold    = number of old nodes
[neltold,~]  = size(triconOld);  % neltold  = number of old elements
[nedgeold,~] = size(edgeconOld); % nedgeold = number of old edges

pNew       = zeros(3, nnold+nedgeold);  % add one node per old edge
triconNew  = zeros(4*neltold, 3);       % four new elts for each old one
edgeconNew = zeros(2*nedgeold + 3*neltold, 2); % each old edge makes 2 new ones,
                                               % plus each elt. makes 3
triedgeNew = zeros(4*neltold, 3);


% Copy old nodes into new list
pNew(:,1:nnold) = pOld;

% Put a new node in the center of each old edge,
% and make the old edge into two new ones
for i = 1:nedgeold
   pn = mean(pOld(:,edgeconOld(i,:)),2);
   if sphere
       pNew(:,nnold+i) = pn/norm(pn);  % Add new node at unit length
   else
       pNew(:,nnold+i) = pn;   % New node directly between old nodes
   end
   edgeconNew(i,:)          = [edgeconOld(i,1), nnold+i];
   edgeconNew(nedgeold+i,:) = [nnold+i, edgeconOld(i,2)];
end

% Generate new elements, and new edges interior to old elements
ne = zeros(9,1);  % storage for new-edge numbers within one element
for i = 1:neltold
    on = triconOld(i,:);  % old (corner) node numbers
    oe = triedgeOld(i,:); % old edge numbers
    nn = nnold+oe;   % new (mid-edge) node numbers
    % Create 3 new edges, interior to the old element
    ne(7:9) = 2*nedgeold + (3*i)-[2,1,0]; % new edge numbers
    edgeconNew(ne(7),:) = [nn(1);nn(2)];
    edgeconNew(ne(8),:) = [nn(2);nn(3)];
    edgeconNew(ne(9),:) = [nn(3);nn(1)];
    % New elements
    triconNew(4*i-3,:) = [on(1); nn(1); nn(3)];
    triconNew(4*i-2,:) = [on(2); nn(2); nn(1)];
    triconNew(4*i-1,:) = [on(3); nn(3); nn(2)];
    triconNew(4*i,  :) = [nn(1); nn(2); nn(3)];
    % Fill in numbers of new edges along edged of old element.
    % Note that the old edges may run in either direction
    for j = 1:3
        if edgeconOld(oe(j),1) == on(j)
            ne(2*j-1) = oe(j);
            ne(2*j)   = oe(j)+nedgeold;
        else
        % Edge connectivity for new elements
            ne(2*j-1) = oe(j)+nedgeold;
            ne(2*j)   = oe(j);
        end
    end
    % Using the new-edge numbers, fill in edge connectivity for new elts.
    triedgeNew(4*i-3,:) = [ne(1); ne(9); ne(6)];
    triedgeNew(4*i-2,:) = [ne(3); ne(7); ne(2)];
    triedgeNew(4*i-1,:) = [ne(5); ne(8); ne(4)];
    triedgeNew(4*i,  :) = [ne(7); ne(8); ne(9)];
end

return;