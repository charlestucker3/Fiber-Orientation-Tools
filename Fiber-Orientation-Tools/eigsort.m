function [lambda, R] = eigsort(A)
% EIGSORT   Eigenvalues and eigenvectors, sorted in descending order
%    [LAMBDA, R] = EIGSORT(A) returns the eigenvalues and eigenvectors of a
%    square matrix A, in descending order. LAMBDA is a row vector of the
%    eigenvalues and R is a matrix the same size as A whose columns contain
%    the components of the eigenvectors.  The original matrix can be
%    recovered from A = R * diag(LAMBDA) * R'.

[Evecs, Aprin] = eig(A);    % unsorted eigenvectors and eigenvalues
% sort eigenvalues in descending order
[lambda, sortOrder] = sort(diag(Aprin), 'descend'); 
% sorted eigenvectors form the rotation matrix
R = Evecs(:, sortOrder);

return