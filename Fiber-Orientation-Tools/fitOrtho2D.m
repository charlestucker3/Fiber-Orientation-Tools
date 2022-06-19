function [a, varargout] = fitOrtho2D(type, np, varargin)
%A = FITORTHO2D(TYPE, NP) generates the polynomial
%    coefficients A such that POLYVAL(A, LAMBDA1) approximates
%    eta1(lambda1) for a planar orthotropic fourth-order orientation
%    tensor.  The resulting coefficients can be used to form a closure
%    approximation. TYPE is a single character that determines the function
%    to be fitted. Choices are:
%
%       B  The Bingham distribution
%       D  Transient distribution function calculation in planar
%          elongation.  In this case a fourth input argument must be
%          present giving the interaction coefficient CI.
%       E  The planar ellipse-radius distribution.
%
%    NP is the order of the polynomial.
%
%A = FITORTHO2D(TYPE, NP, IPLOT) with IPLOT = 1 also plots the data and 
%    fit.  The plot is made in the current figure or subplot. For any other
%    value of IPLOT no plot is created.
%
%    All fits are constrained to eta1 = 0 at lambda1 = 0.5, while the
%    Bingham and ellipse-radius fits are also constrained to eta1 = 0 at
%    lambda1 = 1.
%
%[A, ERR] = FITORTHO2D(TYPE, NP) also returns the RMS error for the data
%    points.
%
%    Fits with NP = 5 are useful for most cases.  If you receive a warning
%    that the matrix is close to singular or badly scaled, try reducing NP.
%
%    See also: CLOSEA4PLANAR, PRINTCOEFF, FITNONORTHO2D

%% Generate the data to be fitted
switch upper(type)
    case 'B'  % Bingham distribution function
        % -- Compute the data to be fitted
        lambda1 = linspace(0.5, 1, 101)'; % Make lambda1 a column vector
        eta1    = zeros(size(lambda1));
        npoints = 180;
        for i = 1:length(lambda1)-1 % Skip computing the last point
            A = diag([lambda1(i), 1-lambda1(i)]);
            [~, ~, ~, ~, A4bing] = fitBingham2D(A, npoints);
            eta1(i) = A4bing(1,1) - 0.5 * (lambda1(i) + lambda1(i)^2);
        end
        ncon = 2;  % Number of constraints for the fit
        
    case 'D' % Distribution function solution, planar elongation
        if nargin >= 4
            CI = varargin{2};
        else
            error('Include CI in input arguments for type == D')
        end
        xi = 1;
        L = [1, 0; 0,-1];  % Velocity gradient: planar elongation
        nnod = 180;        % Points in phi for distribution function calculation
        tend = 5;          % Max simulation time
        dt   = 0.05;       % Simulation time step
        Co   = 0.25;       % Courant number in distribution function calc.
        % Transient distribution function history
        [t, psi, phi] = solvePsi2D(L, CI, xi, tend, dt, nnod, Co);
        p = [cos(phi); sin(phi)];         % p vectors for the set of phi values
        lambda1 = zeros(length(t), 1);    % Make this a column vector
        eta1    = zeros(size(lambda1));
        for i = 1:length(t)
            [A, A4] = p2A(p, psi(:,i));
            lambda1(i) = A(1,1);
            eta1(i)    = A4(1,1) - 0.5*(lambda1(i) + lambda1(i)^2);
        end
        ncon = 1;   % Use only one constraint
        
    case 'E'  % Ellipse radius distribution function
        % -- Compute the data to be fitted
        lambda1 = linspace(0.5, 1, 101)'; % Make lambda1 a column vector
        eta1    = zeros(size(lambda1));
        npoints = 180;
        for i = 1:length(lambda1)-1 % Skip computing the last point
            A = diag([lambda1(i), 1-lambda1(i)]);
            [~, ~, ~, ~, A4er] = fitERdistn2D(A, npoints);
            eta1(i) = A4er(1,1) - 0.5 * (lambda1(i) + lambda1(i)^2);
        end
        ncon = 2;  % Number of constraints for the fit
        
    otherwise
        error('%s is not a legal value of TYPE', type)
end


%%  Fit a polynomial to the data, with constraints
%   Here we follow the notation of Gavin, 
%   https://people.duke.edu/~hpgavin/cee201/constrained-least-squares.pdf,
%   accessed 4/8/2022.
n = np+1;            % Number of fitting parameters. 
T = lambda1.^(np:-1:0);
% The constraints are expressed as b = A*a
if ncon == 1
    % Constrain eta1 = 0 at lambda1 = 0.5 only (for dist'n. fcn. fits)
    b = 0;
    A = zeros(1, n);
    A(:,n) = 1;
    Z = zeros(1,1);
    for i = 1:n-1
        A(:,i) = 0.5.^(n-i);
    end
elseif ncon == 2
    % Constraints are eta1 = 0 at lambda1 = 0.5 and 1 (for Bingham)
    b = zeros(2,1);
    A = zeros(2, n);
    A(:,n) = 1;
    for i = 1:n-1
        A(:,i) = [0.5; 1].^(n-i);
    end
    Z = zeros(2,2);
else
    error('Should never get here')
end

% Solve the constrained least squares problem
x = [T'*T, A'; A, Z] \ [T'*eta1; b];
a = x(1:n);  % Polynomial coefficients; the remaining coefficients
             % in x are the Lagrange multipliers

if nargout >= 2
    % Compute RMS errors at data points
    err = sqrt(mean((eta1 - T*a).^2));
    varargout{1} = err;
end

if nargin >= 3 && varargin{1} == 1
    % Plot data and add lines for the fits
%     figure; clf; 
    hold on
    % Plot the data used for the fit
    plot(lambda1, eta1, 'bs', 'MarkerSize', 8)
    % Plot the fitted polynomial
    lambdaPlot = linspace(0.5, 1, 101);
    plot(lambdaPlot, polyval(a,   lambdaPlot), 'b-')
    % Add the upper and lower bounds on eta1
    plot(lambdaPlot, (-lambdaPlot.^2 + lambdaPlot)/2, 'k-')
    plot(lambdaPlot, ( lambdaPlot.^2 - lambdaPlot)/2, 'k--')
    
    set(gca, 'FontSize', 16)
    xlabel('\lambda_1')
    ylabel('\eta_1')
    legend('Data', 'Fit', 'Upper bound', 'Lower bound' ,...
           'Location', 'NorthWest')
    grid on
    box on
end

end
