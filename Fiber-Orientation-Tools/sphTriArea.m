function [area, varargout] = sphTriArea(pnod, tricon)
%AREA = SPHTRIAREA(PNOD, TRICON) finds the area of each spherical triangle
%    in a triangular mesh.  PNOD (3xnumnod) gives the Cartesian coordinates
%    of each node in the mesh and TRICON (numelt x 3) lists the nodes
%    in each element.  AREA(I) is the spherical area of element I, in the 
%    same order as TRICON.  The P vectors must have unit length.  
%
%[AREA, PCEN] = SPHTRIAREA(PNOD, TRICON) also returns the coordinates PCEN
%    (3xnumelt) of the centroid of each element.

[numelt, ~] = size(tricon);  % number of elements
area = zeros(numelt, 1);     % Storage for element areas
for j = 1:numelt
    % -- Find area of spherical triangle corresponding to element
    %    See Wikipedia article on Spherical Trigonometry
    a = acos(pnod(:,tricon(j,2))' * pnod(:,tricon(j,3)));
    b = acos(pnod(:,tricon(j,3))' * pnod(:,tricon(j,1)));
    c = acos(pnod(:,tricon(j,1))' * pnod(:,tricon(j,2)));
    s = (a+b+c)/2;
    area(j) = 4 * atan(sqrt(tan(s/2) * ...
                            tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)));
end

if nargout >= 2  % Then also calculate and return the element centroids
    pcen = zeros(3, numelt);  %  Storage for element-centroid p values
    for j = 1:numelt
        pcen(:,j) = mean(pnod(:, tricon(j,:)), 2);
        pcen(:,j) = pcen(:,j)/norm(pcen(:,j));
    end
    varargout{1} = pcen;
end

return