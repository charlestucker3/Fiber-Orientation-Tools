function [edgecon, triedge] = meshcon(tricon)
%[EDGECON, TRIEDGE] = MESHCON(TRICON) builds the edge connectivity 
%     information for a triangular mesh, as needed by REFINEMESH.
%     TRICON (numelt x 3) is the element connectivity matrix.  
%     EDGECON (numedge x 2) gives the node numbers for each element edge.
%     TRIEDGE (numelt x 3) gives the edge numbers for each element, in 
%     the same order as the node numbers in TRICON.  TRIEDGE(I,1) 
%     is the number of the edge connecting local nodes 1 and 2 in 
%     element I, TRIEDGE(I,2) connects local nodes 2 and 3, and
%     TRIEDGE(I,3) connects local nodes 3 and 1. 
%
%     See also: REFINEMESH

%% -- Build edgecon
% edgecon(i,j) = local node j in edge i
% THe or rder of nodes for each edge does not matter
[numelt, ~] = size(tricon);  % Number of elements
numnod = max(max(tricon));   % Largest node number

% -- Figure out which nodes are connected via edges
%    isedge(i,j) will equal 1 if nodes i and j share an edge.
isedge = zeros(size(numnod, numnod)); 
for i = 1:numelt
    for j = 1:3
        k = mod(j, 3)+1;  % local number of the second node in the edge
        isedge(tricon(i,j), tricon(i,k)) = 1;
    end
end
numedge = sum(sum(isedge));  % The number of edges in the mesh

% -- Store the edge info in edgecon
edgecon = zeros(numedge, 2);
nextedge = 1;
for i = 1:numnod
    for j = 1:numnod
        if isedge(i,j)
            edgecon(nextedge, 1) = i;
            edgecon(nextedge, 2) = j;
            nextedge = nextedge + 1;
        end
    end
end

%% -- Build triedge
% triedge(i,j) = edge number j in element i.  
% Edge 1 connects nodes tricon(i,1) and tricon(i,2),
% edge 2 connects nodes tricon(i,2) and tricon(i,3), and
% edge 3 connects nodes tricon(i,3) and tricon(i,1).
triedge = zeros(size(tricon));
for i = 1:numelt
    for j = 1:3  
        k = mod(j, 3)+1;  % local number of the second node in the edge
        % Find the row number in edgecon for edge j of element i
        triedge(i,j) = find((edgecon(:,1) == tricon(i,j) & ...
                             edgecon(:,2) == tricon(i,k)) ...
                             | ...
                            (edgecon(:,2) == tricon(i,j) & ...
                             edgecon(:,1) == tricon(i,k)),    1);
    end
end

return
