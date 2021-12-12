function fill3elt(p,tricon,fval)%FILL3ELT(P, TRICON, FVAL) creates a "fill" plot of the element values%    of a scalar function on a 3-D mesh of triangles.  %    P (3 x numnode) gives the nodal coordinates, TRICON (numelt x 3)%    is the nodal connectivity for each element, and FVAL (numelt) is%    a vector of the function values for each element.  %%    See FILL3MESH for a plot using nodal values.hold on% Figure out number of elements & nodes per element[nelts, ~] = size(tricon);% Plot the elementsfor i = 1:nelts    pn = p(:,tricon(i,:)); % Each column of pn gives the x,y,z coords of a node    fill3(pn(1,:), pn(2,:), pn(3,:), fval(i));end% Dress up the plotaxis equalview(3)set(gca, 'FontSize', 16)xlabel('x'); ylabel('y'); zlabel('z');box ongrid on%axis off