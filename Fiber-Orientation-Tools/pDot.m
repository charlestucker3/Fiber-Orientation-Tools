function dpdt = pDot(p, L, xi)
%DPDT = PDOT(P, L, XI) gives the time rate of change of the orientation 
%    vector P (3x1) according to Jeffery's equation, for a velocity
%    gradient tensor L (3x3) and a scalar particle shape factor XI.

W = (L-L')/2;    % Vorticity tensor
D = (L+L')/2;    % Rate-of-deformation tensor
dpdt = W*p + xi*(D*p - (p'*D*p)*p);  % Jeffery's equation

end

