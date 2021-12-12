function [t,p, varargout] = solveJeffery(pzero, L, xi, tspan, varargin)
%[T, P] = SOLVEJEFFERY(PZERO, L, XI, TSPAN) solves Jeffery's equation 
%   for the motion of a single fiber.  PZERO (3x1) is the initial 
%   orientation vector, L (3x3) is the velocity gradient, XI is the 
%   particle shape factor.  The function returns a vector of times T
%   and an array of orientation vectors at those times P.
%
%[T, P, THETA, PHI] = SOLVEJEFFERY(PZERO, L, XI, TSPAN) also returns
%   the polar-coordinate angles THETA and PHI of the vector P at each 
%   time, with 0 <= THETA <= pi and -pi/2 <= PHI <= pi/2.
%
%[T, P] = SOLVEJEFFERY(PZERO, L, XI, TSPAN, OPTIONS) passes the structure
%   OPTIONS (created using the ODESET function) to ODE45.  Otherwise, the 
%   default solution tolerances are RelTol = 1e-6 and AbsTol = 1e-9.  

if nargin == 4
    % Use default solution tolerances
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Options for ode45
else
    % Use user-supplied tolerances
    options = varargin{1};
end

% Primary calculation
[t, p] = ode45(@(t,p) pDot(p, L, xi), tspan, pzero, options);

if nargout >= 3
    % Recover (theta, phi) description
    theta = acos(p(:,3));          % Returns    0  <= theta <= pi
    phi   = atan(p(:,2)./p(:,1));  % Returns -pi/2 <=  phi  <= pi/2
    varargout{1} = theta;
    varargout{2} = phi;
end
