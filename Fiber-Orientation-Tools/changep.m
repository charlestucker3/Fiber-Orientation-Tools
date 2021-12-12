function p = changep(F, pzero)
%p = CHANGEP(F, pzero) Transform a set of initial orientation vectors 
%    pzero (3xn) using pseudo-affine deformation according to a 
%    deformation tensor F (3x3):  p = F dot pzero / |F dot pzero|.
%    Normally we expect det(F) = 1, to create a constant-volume deformation.
%    Returns p, which is the same size as pzero.

p = F*pzero;  % Affine deformation; new p's don't have unit length

% Lengths of the p vectors from above
pnorm = sqrt(sum(p.*p));
% Normalize the lengths of the new vectors
p(1,:) = p(1,:)./pnorm;
p(2,:) = p(2,:)./pnorm;
p(3,:) = p(3,:)./pnorm;

return;