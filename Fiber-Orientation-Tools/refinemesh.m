function [pNew, triconNew, varargout] = ...
             refinemesh(pOld, triconOld, edgeconOld, triedgeOld, nseg)
%[PNEW, TRICONNEW] = REFINEMESH(P, TRICON, EDGECON, TRIEDGE, NSEG) 
%   refines a triangular mesh by dividing each edge into NSEG segments and
%   building a uniform mesh with NSEG^2 triangular elements within each old
%   element.  All new nodes are in the plane of the old element. 
%
%   P (3xn) contains the old nodal coordinates, TRICON (nelt x 3) lists the
%   nodes for each element, EDGECON (nedge x 2) lists the nodes for each
%   edge and TRIEDGE (nelt x 3) lists the edges in each element. The order
%   of edges in TRIEDGE must match the order of nodes in TRICON: Edge 1
%   connects local nodes 1 and 2, edge 2 connects local nodes 2 and 3, and
%   edge 3 connects local nodes 3 and 1.%   EDGECON and TRIEDGE can be
%   generated from TRICON using the MESHCON function.
%
%   PNEW and TRICONNEW give the nodal coordinates and connectivity of the
%   new mesh.
%
%[PNEW, TRICONNEW, PCEN] = REFINEMESH(... also returns PCEN (3 x newnodes)
%   giving the element centroid coordinates.
%
%   See also: MESHCON


if nseg < 1
    error('NSEG must be >= 1')
end

if nseg == 1 % Special case: no mesh refinement
    pNew = pOld;
    triconNew = triconOld;
    if nargout >= 3
        % Compute and return element centroids
        [neltold,~]  = size(triconOld);
        pcen = zeros(3, neltold);
        for i = 1:neltold
            pcen(:,i) = mean(pNew(:,triconNew(i,:)), 2);
        end
        varargout{1} = pcen;
    end
    return
end

%% Old sizes
[~,nnold]    = size(pOld);       % nnold    = number of old nodes
[neltold,~]  = size(triconOld);  % neltold  = number of old elements
[nedgeold,~] = size(edgeconOld); % nedgeold = number of old edges

% New sizes and storage space
nnnew      = nnold + (nseg-1)*nedgeold + (nseg-1)*(nseg-2)*neltold/2; 
%          = number of new nodes
neltnew    = neltold * nseg^2;      % number of new elements
pNew       = zeros(3, nnnew); 
triconNew  = zeros(neltnew, 3); 
childnode  = zeros(nedgeold, nseg+1); 
%          = the node numbers in the new mesh along each old edge
edgevec    = zeros(3, nedgeold); % Unit vector along each old edge
dl         = zeros(nedgeold,1);  % New elt. edge length for each old edge 


% Copy old nodes into new list
pNew(:,1:nnold) = pOld;

% Add (nseg-1) new nodes to each old edge and compute unit vectors for edges 
for i = 1:nedgeold
    % Edge length and unit vector from first edge node to second
    edgevec(:,i) = pOld(:,edgeconOld(i,2)) - pOld(:,edgeconOld(i,1));
    edegelength = norm(edgevec(:,i));         % Edge length
    edgevec(:,i) = edgevec(:,i)/edegelength;  % Make this a unit vector
    dl(i) = edegelength/nseg;     % Distance between new nodes on this edge
    
    % New node numbers and nodal coordinates
    newnum = nnold + (i-1)*(nseg-1) + (1:nseg-1);  % Vector of new numbers
    pNew(:, newnum) = pOld(:,edgeconOld(i,1)) ...
                    + edgevec(:,i)*(1:nseg-1)*dl(i);
    
    % Child nodes for this edge include the two old end nodes
    % and the new nodes just created
    childnode(i,:) = [edgeconOld(i,1), newnum, edgeconOld(i,2)];
end

% Generate new nodes and elements within each old element
for i = 1:neltold
    % Preliminary quantities for this element
    % - List of nodes along "left" edge, from node 1 to node 2
    leftedge = childnode(triedgeOld(i, 1),:);
    if leftedge(1) ~= triconOld(i,1)
        leftedge = fliplr(leftedge); 
        % Flip the list, if necessary to get the proper order
    end
    % - List of nodes along "right" edge, from node 1 to node 3
    rightedge = childnode(triedgeOld(i, 3),:);
    if rightedge(1) ~= triconOld(i,1)
        rightedge = fliplr(rightedge); 
    end
    % - List of nodes across "bottom" edge, from node 2 to node 3
    bottomedge = childnode(triedgeOld(i,2),:);
    dp = edgevec(:,triedgeOld(i,2));  % Vector along this edge
    if bottomedge(1) ~= triconOld(i,2)
        bottomedge = fliplr(bottomedge); 
        dp = -dp;  % Also flip the dp vector in this case
    end
    % Set length of dp to the edge length of new nodes in this direction
    dp = dp * dl(triedgeOld(i,2));
    
    % Generate new internal nodes and new elements in rows, starting
    % from node 1 and working down.
    prevrow  = zeros(1, nseg-1); % Node numbers at bottom of previous row
                                 % = node numbers at top of this row
    newrow   = zeros(1, nseg);   % Node numbers at bottom of this row
    lastnod = nnold + (nseg-1)*nedgeold + (nseg-1)*(nseg-2)*(i-1)/2;
    %       = index of last filled spot in pNew
    lastelt = (i-1)*nseg^2;  % Index of last filled spot in triconNew
    for j = 1:nseg
        if j ==1
            % Top row has one element and no new nodes
            prevrow(1)  = leftedge(1);
            newrow(1:2) = [leftedge(2), rightedge(2)];
            triconNew(lastelt+1,:) = [prevrow(1), newrow(1), newrow(2)];
            lastelt = lastelt + 1;
        else
            % Interior rows
            % First set up the lists of nodes for this row
            prevrow = newrow;  % The top line of nodes for this row
            if j < nseg
                % Create new interior nodes for the lower lines of nodes
                newrow = [leftedge(j+1), lastnod + (1:j-1), rightedge(j+1)];
                pNew(:, lastnod + (1:j-1)) = pNew(:,newrow(1)) ...
                                           + dp * (1:j-1);
                lastnod = lastnod + j-1;
            else
                % For the last row, use the nodes along the bottom edge
                % of the parent element
                newrow = bottomedge;
            end
            
            % Then generate two lines of elements along this row
            for k = 1:j-1
                % Elements in the first line "hang" from the top nodes
                triconNew(lastelt + k, :) = ...
                             [prevrow(k); newrow(k+1); prevrow(k+1)];
            end
            lastelt = lastelt + j-1;
            
            for k = 1:j
                % Elements in the second line "sit" on the bottom nodes
                triconNew(lastelt + k, :) = ...
                             [newrow(k); newrow(k+1); prevrow(k)];
            end
            lastelt = lastelt + j;
        end
    end
        
end

% Calculated element centroid coordinates, if requested
if nargout >= 3
    pcen = zeros(3, neltnew);
    for i = 1:neltnew
        pcen(:,i) = mean(pNew(:,triconNew(i,:)), 2);
    end
    varargout{1} = pcen;
end

return;