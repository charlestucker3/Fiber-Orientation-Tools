function [a, varargout] = polyfit2(x, y, z, np, varargin)
%A = POLYFIT2(X, Y, Z, NP) will fit a full 2-D polynomial of order NP, 
%    with independent variables X and Y, to the points Z.  X, Y and Z can
%    be vectors or matrices, but must be the same size.  The coefficients A
%    are returned in the order used by POLYVAL2.
%
%A = POLYFIT2(X, Y, Z, NP, XC, YC, ZC) will constrain the polynomial to 
%    pass through values ZC at points (XC, YC).  XC, YC and ZC must be
%    vectors of the same length.
%
%[A ERR] = POLYFIT2(X, Y, Z, NP) also returns the RMS error at the data
%    points.  
%
%    Currently NP <= 5, though this could be extended. 
%
%    See also: POLYVAL2

% Ensure that the input arrays are the same size
[xrow, xcol] = size(x);
[yrow, ycol] = size(y);
[zrow, zcol] = size(z);
if yrow ~= xrow || ycol ~= xcol || zrow ~= xrow || zcol ~= xcol
    error('X, Y and Z must be the same size and shape')
end

% Make x, y and z into column vectors
ndat = xrow * xcol;       % Number of data points
x = reshape(x, ndat, 1);
y = reshape(y, ndat, 1);
z = reshape(z, ndat, 1);


if np > 5
    error('Currently NP must be <= 5')
end

% The l and m vectors are the powers of x and y associated with each term.
% To treat higher-order polynomials, extend l and m
l = [0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0];
m = [0 0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5];
% Trim l and m to the length needed for a polynomical of order np
n = (np+1) * (np+2) / 2;  % Number of active coefficients
l = l(1:n);
m = m(1:n);

% Form the matrix T, such that z = T * a gives fitted values at the
% data points (x, y).
 T = x.^l .* y.^m;  

if nargin == 4
    % Do an unconstrained fit
    a = (T'*T) \ (T'*z);
    
elseif nargin >= 7
    % Do a constrained fit
    % Store and check the input data
    xc = varargin{1};
    yc = varargin{2};
    zc = varargin{3};
    nc = length(xc);
    if length(yc) ~= nc || length(zc) ~= nc
        error('XC, YC and ZC must have the same length')
    end
    % Ensure that xc, yc and zc are column vectors
    xc = reshape(xc, nc, 1);
    yc = reshape(yc, nc, 1);
    zc = reshape(zc, nc, 1);
    % Form the matrixs A, such that the constraints are satisfied 
    % if zc = A * a.
    A = xc.^l .* yc.^m;
    Z = zeros(nc, nc);
    % Compute the constrained fit
    b = [T'*T, A'; A, Z] \ [T'*z; zc];
    a = b(1:n);  % Polynomial coefficients; the remaining coefficients 
                 % in b are the Lagrange multipliers
        
else
    error('Input arguments for constraints are incomplete')
end
        
if nargout >= 1
    % Return the RMS error at the data points
    zfit = T * a;
    err = sqrt(mean((zfit - z).^2));
    varargout{1} = err;
end

end