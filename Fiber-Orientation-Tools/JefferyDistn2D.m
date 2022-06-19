function [psi, phi, varargout] = JefferyDistn2D(A, npoints)
%[PSI, PHI] = JEFFERYDISTN2D(A, NPOINTS) returns the planar Jeffery 
%    orientation distribution function PSI corresponding to the
%    second-order orientation tensor A (2x2).  NPOINTS is the number of
%    equally-spaced points in the interval  -pi/2 <= PHI < pi/2 where the
%    distribution function is evaluated.
%
%[PSI, PHI, A, A4] = JEFFERYDISTN2D(A, NPOINTS) also returns the 
%    second-order tensor A (2x2) and fourth-order tensor A4 (5x5)
%    corresponding to the values PSI and PHI.  The A tensor on output may
%    differ slightly from the input A tensor, owing to the discrete nature
%    of PSI.  These errors can become large for lambda1 > 0.9 unless
%    npoints is very large.  The original value of A and A4 returned by
%    closeA4planar.m are much more accurate.


phi  = linspace(-pi/2, pi/2, npoints+1);
phi  = phi(1:end-1);           % phi values, in [-pi/2, pi/2)
dphi = pi/npoints;             % Delta(phi)

pp = [cos(phi).*cos(phi); ...  % Each column gives the components of the
      sin(phi).*sin(phi); ...  % pp tensor for one phi value.
      sin(phi).*cos(phi)]; 
  
[Eval, Q] = eigsort(A);     % Eigenvalues and eigenvectors of A, sorted


% Form B-inverse in the principal coordinates of A
Binv = diag([ (Eval(2)/Eval(1)), (Eval(1)/Eval(2))] );
%  Then rotate it to the laboratory axes
Binv = Q * Binv * Q';

% Planar Jeffery distribution
psi = 1 ./ (2*pi * (tens2vec(Binv)' * diag([1,1,2]) * pp));


if nargout >= 3
    if Eval(1) > 0.97
        % Accuracy also depends on npoints
        fprintf('WARNING: JefferyDistn2D tensors may not be accurate\n')
    end
    Av = 2 * pp * psi' * dphi;  % Contracted form
    A  = vec2tens(Av);          % Matrix form
    varargout{1} = A;
end

if nargout >= 4
    % Calculate the fourth-order tensor
    p4 = [cos(phi).^4; ...     % Each column of p4 contains the components
          sin(phi).^4; ...     % of the pppp tensor for one phi value.
          cos(phi).^2 .* sin(phi).^2; ...
          cos(phi).^3 .* sin(phi); ...
          cos(phi)    .* sin(phi).^3];

      A4v = 2 * p4 * psi' * dphi;     % Fourth-order tensor, contracted
      A4  = vec2tens4(A4v);           % Matrix form
      varargout{2} = A4;
end

end

