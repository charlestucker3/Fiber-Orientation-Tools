function z = polyval2(a, x, y)
%Z = POLYVAL2(A, X, Y) returns the value Z of a 2-D polynomial in X and Y
%    with coefficients A, as created by POLYFIT2.  X and Y can be scalars
%    or they can be matrices of the same size.  Z with be the same size as
%    X and Y.

n = length(a);  % The number of terms in the polynomial
if n > 21
    error('polyval2 is currently limited to NP <= 5')
end

% If x and y are matrices, reshape them into column vectors.
[rows, cols] = size(x);
[yr,   yc]   = size(y);
if yr ~= rows || yc ~= cols
    error('X and Y must be the same size')
end
x = reshape(x, rows*cols, 1);
y = reshape(y, rows*cols, 1);

% a must also be a column vector
if isrow(a)
    a = a';
end

% The l and m vectors are the powers of x and y associated with each term.
% To treat higher-order polynomials, extend l and m
l = [0 1 0 2 1 0 3 2 1 0 4 3 2 1 0 5 4 3 2 1 0];
m = [0 0 1 0 1 2 0 1 2 3 0 1 2 3 4 0 1 2 3 4 5];
% Trim l and m to the length needed for the given polynomical
l = l(1:n);
m = m(1:n);

% Compute z, initially as a column vector
z = (x.^l .* y.^m) * a;

% Reshape z to match the input shapes of x and y
z = reshape(z, rows, cols);

end
