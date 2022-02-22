function p = changep(F, pprime)
%P = CHANGEP(F, PPRIME) Transforms a set of initial orientation vectors 
%    PPRIME (3xn) according to a deformation tensor F (3x3), using the
%    deformation form of Jeffery's equation for shape factor = 1.
%    Normally det(F) = 1, to create a constant-volume deformation.
%    Returns the final orientation P, which is the same size as PPRIME.

p = F * pprime;  % Affine deformation; new p's don't have unit length

% Lengths of the p vectors from above
pnorm = vecnorm(p);
% Normalize the lengths of the new vectors
p(1,:) = p(1,:)./pnorm;
p(2,:) = p(2,:)./pnorm;
p(3,:) = p(3,:)./pnorm;

return;