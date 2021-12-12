function [A, varargout] = p2A(p, varargin)
%P2A Convert orientation vectors to an orientation tensor.
%   A = P2A(P) gives the second-order orientation tensor A, 
%   in 3x3 matrix form, for the set of orientation vectors P.
%   P is 3xN, with P(I,J) giving the Ith Cartesian component of vector J.
%
%   A = P2A(P, W) weights each vector P(:,J) with the weighting factor
%   W(J). 
%
%[A, A4] = P2A(P) also returns the fourth-order orientation tensor A4,
%   in 6x6 matrix form.

% -- Set up the weighting factors w as a row vector
[~, nump] = size(p);  % The number of p vectors (also used for A4)
if nargin > 1
    % Use the second argument as the weighting factors
    w = varargin{1};
    [wrows, ~] = size(w);
    if wrows > 1, w = w'; end   % Make sure w is a row vector
else
    w = ones(1, nump);  % The default is equal weithing
end

% -- Second-order tensor
A = zeros(3);
% Form the weighted sums
A(1,1) = sum(p(1,:).*p(1,:).*w) / sum(w);
A(2,2) = sum(p(2,:).*p(2,:).*w) / sum(w);
A(3,3) = sum(p(3,:).*p(3,:).*w) / sum(w);
A(2,3) = sum(p(2,:).*p(3,:).*w) / sum(w);
A(3,1) = sum(p(3,:).*p(1,:).*w) / sum(w);
A(1,2) = sum(p(1,:).*p(2,:).*w) / sum(w);
% Fill in remaining components by symmetry
A(3,2) = A(2,3);
A(1,3) = A(3,1);
A(2,1) = A(1,2);

% -- Fourth-order tensor
if nargout >= 2  
    A4 = zeros(6,6);      % Storage for the fourth-order tensor
    wtot = sum(w);        % Saves re-computing the sum
    for i = 1:nump
        pp  = p(:,i) * p(:,i)';  % The dyadic produce of p with itself
        ppv = tens2vec(pp);      % pp in 6x1 column-vector form
        A4 = A4 + w(i)*(ppv*ppv') / wtot; % ppv*ppv' = pppp, 6x1 form
    end
    varargout{1} = A4;
end

return