function [A, B, D, varargout] = Clayer2laminate(Clayer, z)
%[A, B, D] = Clayer2laminate(Clayer, Z) computes the laminate
%    stiffness matrices A, B and D (all 3x3) for a N-layer laminate where
%    the 3-D stiffness tensor of the Ith layer is contained in
%    Clayer(:,:,I) and layer I runs from Z(I) to Z(I+1), with Z=0 on the
%    midplane. The coordinate system for CLAYER must have Z as the
%    thickness direction of the laminate.
%
%    Classical Kirchoff-Love laminated plate theory is used.
%
%[A, B, D, ETENSILE, EFLEX] = Clayer2laminate(Clayer, Z) also returns the
%    tensile and flexural moduli, with ETENSILE(1) and EFLEX(1) being the
%    X-direction properties and ETENSILE(2) and EFLEX(2)the properties in
%    the Y direction.
%

% Total laminate thickness
H = z(end) - z(1);

% Check that z = 0 on the midplane
tiny = 1e-9;
if abs(z(1)+H/2) > tiny || abs(z(end)-H/2) > tiny
    error('z values must range from -H/2 to H/2')
end

% Number of layers
[~, ~, n] = size(Clayer);
if length(z) ~= n+1
    error('length(z) must be number of layers + 1')
end

% Storage for final results
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);
if nargout >= 4
    Etensile = zeros(2,1);
    Eflex    = zeros(2,1);
end

% Loop over the layers and build up the A, B and D matrices
for i = 1:n
    S = inv(Clayer(:,:,i));  % Traditional compliance matrix for layer i
    Q = inv([S(1,1), S(1,2), S(1,6); ...
             S(2,1), S(2,2), S(2,6); ...
             S(6,1), S(6,2), S(6,6)]);  % Plane-stress stiffness, layer i
    A = A + Q * (z(i+1)   - z(i));
    B = B + Q * (z(i+1)^2 - z(i)^2)/2;
    D = D + Q * (z(i+1)^3 - z(i)^3)/3;
end
% A, B and D are now complete.  

if nargout >= 4  % Compute the tensile and flexural moduli
    Slam = inv([A, B; ...
                B, D]);    % Full compliance matrix for laminate
    Etensile(1) =  1 / (Slam(1,1)*H);
    Etensile(2) =  1 / (Slam(2,2)*H);
    Eflex(1)    = 12 / (Slam(4,4)*H^3);
    Eflex(2)    = 12 / (Slam(5,5)*H^3);
    varargout{1} = Etensile;
    if nargout >= 5
        varargout{2} = Eflex;
    end
end