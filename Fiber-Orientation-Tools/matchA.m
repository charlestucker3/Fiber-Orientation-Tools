function [f, varargout] = matchA(fhat, p, A)
%F = MATCHA(FOLD, P, A), where A is a 3x3 second-order orientation tensor,
%    finds the weights F associated with the orientation vectors P to match
%    the tensor A, while minimizing the mean square deviation between F and
%    FHAT.  This can be used to fine-tune an initial discrete
%    reconstruction of a distribution function that generates FOLD and P.
%    P is 3xn, and both FOLD and F are 1xn.  
%
%[F, NNEG] = MATCHA(FHAT, P, A) also returns the number of negative entries
%    in F as NNEG.  NNEG = 0 indicates that all entries in F are >= 0.  If
%    the NNEG output argument is missing and some of the weights are
%    negative, a warning is printed to the console.  Negative weights can
%    arise when the orientation state has eigenvalues close to zero or
%    unity and the points P are not sufficiently dense.  
%
%F = MATCHA(FHAT, P, A4), where A4 is a 6x6 fourth-order orientation
%    tensor, returns weights F that match A4.  This guarantees that the
%    weights also match the second-order tensor corresponding to A4. 

%    The function uses a minimum-distance solution, minimizing norm(f-fhat)
%    subject to C*f = tensevec(A).  
%    See http://www.seas.ucla.edu/~vandenbe/133A/lectures/cls.pdf, 1/7/2021

% Ensure that fhat is a column vector
fhat = reshape(fhat, length(fhat), 1);

% Determine whether the match will be second-order or fourth-order
[rows, cols] = size(A);

if rows == 3 && cols == 3 % A is a second-order tensor
    % C matrix, such that C*f = tens2vec(A)
    C = [p(1,:).*p(1,:); p(2,:).*p(2,:); p(3,:).*p(3,:); ...
         p(2,:).*p(3,:); p(3,:).*p(1,:); p(1,:).*p(2,:)];
    % QR decomposition of C-transpose
    [Q, R] = qr(C');
    % Final weights
    f = fhat + (Q * (R'\(tens2vec(A) - C*fhat)));
    
elseif rows == 6 && cols == 6 % A is a fourth-order tensor
    % C matrix, such that C*f = tens2vec4(A)
    C = [p(1,:).*p(1,:).*p(1,:).*p(1,:); ...
         p(2,:).*p(2,:).*p(2,:).*p(2,:); ...
         p(3,:).*p(3,:).*p(3,:).*p(3,:); ...
         p(2,:).*p(2,:).*p(3,:).*p(3,:); ...
         p(3,:).*p(3,:).*p(1,:).*p(1,:); ...
         p(1,:).*p(1,:).*p(2,:).*p(2,:); ...
         p(1,:).*p(1,:).*p(2,:).*p(3,:); ...
         p(1,:).*p(1,:).*p(3,:).*p(1,:); ...
         p(1,:).*p(1,:).*p(1,:).*p(2,:); ...
         p(2,:).*p(2,:).*p(2,:).*p(3,:); ...
         p(2,:).*p(2,:).*p(3,:).*p(1,:); ...
         p(2,:).*p(2,:).*p(1,:).*p(2,:); ...
         p(3,:).*p(3,:).*p(2,:).*p(3,:); ...
         p(3,:).*p(3,:).*p(3,:).*p(1,:); ...
         p(3,:).*p(3,:).*p(1,:).*p(2,:)];
    % QR decomposition of C-transpose
    [Q, R] = qr(C');
    % Final weights
    f = fhat + (Q * (R'\(tens2vec4(A) - C*fhat)));
    
else
    error('A must be either 3x3 or 6x6')
end

% Count any negative weights
nneg = sum(f<0);  % nneg is the number of negative weights

if nargout >= 2
    % If there's a second output argument, return nneg
    varargout{1} = nneg;
else
    if nneg > 0
        % If there's no second output argument and there are negative
        % weights, print a warning.
        fprintf('WARNING: %i of %i p vectors have negative weights\n', ...
            nneg, length(p))
    end
end