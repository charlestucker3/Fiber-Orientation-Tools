function A = thetaphi2A(theta, phi, varargin)
%A = THETAPHI2A(THETA, PHI) returns the orientation tensor A (3x3)
%     corresponding to the fiber orientation angles THETA and PHI observed
%     on a planar cross-section.  THETA and PHI are vectors of the same
%     size, and are in radians.  The section plane is assumed to be the 1-2
%     plane. By default, the calculation uses the Bay weighting function
%     with theta_c corresponding to L/d = 20.
%
%A = THETAPHI(THETA, PHI, WTYPE, ASPECT) uses the weighting function WTYPE.
%     If ASPECT is present it is used as the fiber aspect ration L/d in the
%     weighting function.  Otherwise, ASPECT = 20.  Options for WTYPE are
%        'Bay', weighting function is max(1/cos(theta), 1/cos(theta_c))
%                    with cos(theta_c) = d/L
%        'Konicek', weighting function is 1/(cos(theta) + sin(theta)*d/L)
%
%     The components A(1,3) = A(3,1) and A(2,3) = A(3,2) are ambiguous in
%     this measurement.  The function returns the upper bound for these
%     components.  The lower bound has the same magnitude and opposite
%     sign.  

% Parse input arguments
if nargin >= 3
    wtype = varargin{1};
else
    wtype = 'Bay';  % Default weighting function
end
if nargin >= 4
    aspect = varargin{2};
else
    aspect = 20;   % Default fiber aspect ratio for weighting
end

% Calculate the weighting factor for each fiber
switch lower(wtype)
    case 'bay'
        weight = 1./cos(theta);
        weight(weight > aspect) = aspect;
    case 'konicek'
        weight = 1./(cos(theta) + sin(theta)/aspect);
    otherwise
        error('Illegal value of wtype. Try ''Bay'' or ''Konicek''.')
end

 
% Orientation vectors for all fibers (3xn)
p = [cos(phi).*sin(theta); sin(phi).*sin(theta); cos(theta)];

% Build the orientation tensor
A = zeros(3);   
A(1,1) = sum(p(1,:).*p(1,:).*weight) / sum(weight);
A(1,2) = sum(p(1,:).*p(2,:).*weight) / sum(weight);
A(1,3) = sum(abs(p(1,:).*p(3,:)).*weight) / sum(weight);
A(2,1) = A(1,2);
A(2,2) = sum(p(2,:).*p(2,:).*weight) / sum(weight);
A(2,3) = sum(abs(p(2,:).*p(3,:)).*weight) / sum(weight);
A(3,1) = A(1,3);
A(3,2) = A(2,3);
A(3,3) = sum(p(3,:).*p(3,:).*weight) / sum(weight);

return




