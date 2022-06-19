function [a1, a2, err1, err2] = fitNonOrtho2D(CI, xi, n1, n2, ...
                                           alpha, dt, tend, varargin)
%[A1, A2, ERR1, ERR2] = FITNONORTHO(CI, XI, N1, N2, ALPHA, DT, TEND) finds  
%    coefficient vectors A1 and A2 for a planar non-orthotropic closure
%    approximation for flow-induced fiber orientation for the interaction
%    coefficient CI and the particle shape factor XI.
%
%    A1 and A2 are coefficients of complete polynomials in (lambda1, phid)
%    approximating eta1 and eta2.  The polynnomials have orders N1 and N2,
%    respectively. N1 = N2 = 3 is typical.  Setting N1 or N2 too small will
%    generate a warning for a badly scaled or singular matrix.  Use
%    polyval2.m to evaluate the polynomials. ERR1 and ERR2 are the RMS
%    errors in fitting eta1 and eta2.
%
%    ALPHA is a vector describing the flow cases used to generate the
%    fitting data.  For case I, the velocity gradient tensor is L = [0, 1;
%    ALPHA(I), 0].  Thus, ALPHA(I) = 0 is simple shear, ALPHA(I) = 1 is
%    planar elongation, and ALPHA(I) = -1 is rigid-body rotation.  A
%    typical set is ALPHA = [-0.8, -0.15, 0, 0.15, 0.3, 0.5, 0.7, 1];  
%    DT and TEND are the time step and ending time for the fitting data.
%    Values of DT = 0.10 and TEND = 15 are typical, though longer TEND
%    values should be used when CI < 0.05. 
%
%[A1, A2] = FITNONORTHO(CI, XI, N1, N2, ALPHA, DT, TEND, IPLOT) with 
%     IPLOT = 1 generates plots of the fitting data and the fitted
%     functions.  For other values of IPLOT no plots are generated.
%
%[A1, A2] = FITNONORTHO(CI, XI, N1, N2, ALPHA, DT, TEND, IPLOT, CONSTRAINT) 
%     controls constraints for the fitting functions.  If the CONSTRAINT
%     argument is absent, or if CONSTRAINT = 'M', the minimum constraints
%     are applied: eta1 is constrained to be zero at lambda1 = 0.5, phid =
%     0, while eta2 is constrained to be zero along phid = 0 for all values
%     of lambda1. If CONSTRAINT = 'C' then, in addition, eta1 and eta2 = 0
%     for all values of phid along the lambda1 = 0.5 line.
%
%     See also: POLYVAL2, POLYFIT2, FITORTH2D, CLOSEA4PLANAR, PRINTCOEFF

%% Run a set of flow histories and plot their intersections with several
%  planes of constant lambda1.


% Parameters for the distribution function calculations
nnod = 180;        % Points in phi for distribution function calculation
Co   = 0.25;       % Courant number in distribution function calc.

% Parameters controlling the appearance of the polots
% Line colors for each alpha in the fitting data
lcolor = ['c', 'g', 'm', 'r', 'b', 'y', 'k', ...
          'c', 'g', 'm', 'r', 'b', 'y', 'k']; 
fontsize = 16;     % Font size for plot labels, points
linewidth = 2;     % Line width for fitting data
plotsymbol = '.';  % Symbol for fitting data points

%% Generate the fitting data
nalpha = length(alpha);      % Number of alpha values
nt     = round(tend/dt) + 1; % Number of time steps in dist'n. function calc.
% Storage for primary results:
lambda1 = zeros(nalpha, nt);
eta1    = zeros(nalpha, nt);
eta2    = zeros(nalpha, nt);
phid    = zeros(nalpha, nt);

% Loop over the alpha values
for i = 1:nalpha
    L = [0, 1; alpha(i), 0]; % Velocity gradient for this flow case
    [~, dvec] = eigsort(0.5*(L + L'));  % Sorted eigenvectors of D
    d = dvec(:,1);      % Eigenvector of principal eigenvalue of D
    
    % Transient distribution function history
    [t, psi, phi] = solvePsi2D(L, CI, xi, tend, dt, nnod, Co);
    
    % Analyze results for each time step
    p = [cos(phi); sin(phi)];  % p vectors for the phi values
    for j = 1:length(t)
        [A, A4] = p2A(p, psi(:,j));  % Tensors for this time step
        [Evals, Q] = eigsort(A);     % Sorted eigenvalues and eigenvectors
        lambda1(i,j) = Evals(1);     % Largest eigenvalue of A
        A4bar = rotate4(A4, Q');     % A4 in principal axes of A
        eta4 = A4bar - closeA4planar(diag(Evals), 'N');
        eta1(i,j) = eta4(1,1);       % eta1 and eta2 are parameters for ...
        eta2(i,j) = eta4(1,3);       % ... A4bar - A4Natural
        dbar = Q' * d;               % d in principal coords. of A
        phid(i,j) = atan(dbar(2) / dbar(1)); % phi associated with dbar
        % Keep phid within limits, using pi/2 periodicity
%         if phid(i,j) >  pi/4, phid(i,j) = phid(i,j) - pi/2; end
%         if phid(i,j) < -pi/4, phid(i,j) = phid(i,j) + pi/2; end
    end
end
phid(:,1) = 0;  % Patch up phid values for the initial point

%% Fit polynomials in (lambda1, phid) for eta1 and eta2
% -- Determine the constraints to use
if nargin <= 8
    % eta1 is constrained only at the isotropic point
    lambda1c = 0.5;
    phi1c    = 0;
    eta1c    = 0;
    % eta2 is constrained at the isotropic point and along phid = 0 axis
    % The number of points equals n2+1, to avoid duplicate constraints
    lambda2c = linspace(0.5, 1, n2+1); 
    phi2c    = zeros(1, n2+1);
    eta2c    = zeros(1, n2+1);
elseif strcmpi(varargin{2}, 'C')
    % eta1 is constrained at isotropic lambda1 for multiple values of phid
    lambda1c = [0.5, 0.5, 0.5, 0.5];
    phi1c    = [0,   0.5, 1,   1.5];
    eta1c    = [0,   0,   0,   0];
    % eta2 is also constrained at lambda1 = 0.5 for multiple phid values.  
    lambda2c = [0.5, 0.5, 0.5, 0.5, 0.6, 0.8, 1.0];
    phi2c    = [0,   0.5, 1,   1.5, 0,   0,   0  ];
    eta2c    = [0,   0,   0,   0,   0,   0,   0  ];
else
    error('Illegal value of constraint argument')
end
% -- Constrained fits 
[a1, err1] = polyfit2(lambda1, phid, eta1, n1, lambda1c, phi1c, eta1c);
[a2, err2] = polyfit2(lambda1, phid, eta2, n2, lambda2c, phi2c, eta2c);

%% Plotting of data and fitted functions, if desired
if nargin >= 8 && varargin{1} == 1 

    % -- Compute polynomial surfaces for the plots of eta(lambda1, phid)
    [lambdag, phig] = meshgrid(linspace(0.5, 1, 21), linspace(0, pi/2, 21));
    eta1g = polyval2(a1, lambdag, phig);  % Points on the fitted surface
    eta2g = polyval2(a2, lambdag, phig);  % Points on the fitted surface
    % Set values outside the physical range to NaN,
    % so that that portion of the surface will not appear
    [eta1g, eta2g] = limitEta(eta1g, eta2g, lambdag);

    % Plot all the fitting data in lambda1/eta1/eta2 space
    figure(1); clf; hold on
    for i = 1:nalpha
        plot3(lambda1(i,:), eta2(i,:), eta1(i,:), [lcolor(i), '-'], ...
            'LineWidth',  linewidth)
    end
    plot3([0.5,1], [0,0], [0,0], 'k-')  % Central axis
    set(gca, 'FontSize', fontsize)
    xlabel('\lambda_1')
    ylabel('\eta_2')
    zlabel('\eta_1')
    grid on
    box on
    view(60, 30)
%     title('Fitting Data')

    
    % Plot eta1 vs (lambda1, phid) for the fitting data; add fitted surface
    figure(2); clf; hold on
    % Plot fitting data as points
    for i = 1:nalpha
        plot3(lambda1(i,:), phid(i,:), eta1(i,:), [lcolor(i), plotsymbol])
    end
    h1 = surf(lambdag, phig, eta1g, eta1g); % Add the fitted surface
    set(h1, 'FaceAlpha', 0.5)
    colormap pink
%     axval = axis;
%     % Add lines for the limiting values of eta1(lambda1)
%     lambdalimit = linspace(0.5, 1, 51);
%     etalimit    = 0.5 * (lambdalimit - lambdalimit.^2);
%     plot3(lambdalimit, zeros(size(etalimit)),  etalimit, 'k-', ...
%         lambdalimit, zeros(size(etalimit)), -etalimit, 'k-', ...
%         'LineWidth', 1.5)
%     axis(axval)
    set(gca, 'FontSize', fontsize)
    xlabel('\lambda_1')
    ylabel('\phi_d')
    zlabel('\eta_1')
    grid on
    box on
    view(60, 30)
%     title('\eta_1(\lambda_1, \phi_d)')
    
    
    figure(3); clf; hold on
    % Plot fitting data as points
    for i = 1:nalpha
        plot3(lambda1(i,:), phid(i,:), eta2(i,:), [lcolor(i),plotsymbol])
    end
    % Add the fitted surface
    h2 = surf(lambdag, phig, eta2g, eta2g);
    set(h2, 'FaceAlpha', 0.5)
    colormap pink
    set(gca, 'FontSize', fontsize)
    xlabel('\lambda_1')
    ylabel('\phi_d')
    zlabel('\eta_2')
    grid on
    box on
%     title('\eta_2(\lambda_1, \phi_d)')
    view(-125, 25)

end

end



% % % 
% % % figure(6); hold on
% % % h1 = surf(lambdag, phig, eta1g, eta1g);
% % % set(h1, 'FaceAlpha', 0.5)
% % % colormap pink
% % % axval = axis;
% % % % Add lines for the limiting values of eta1(lambda1)
% % % lambdalimit = linspace(0.5, 1, 51);
% % % etalimit    = 0.5 * (lambdalimit - lambdalimit.^2);
% % % plot3(lambdalimit, zeros(size(etalimit)),  etalimit, 'k-', ...
% % %       lambdalimit, zeros(size(etalimit)), -etalimit, 'k-', ...
% % %       'LineWidth', 1.5)
% % % axis(axval)
% % % 
% % % 
% % % figure(7); hold on
% % % h2 = surf(lambdag, phig, eta2g, eta2g);
% % % set(h2, 'FaceAlpha', 0.5)
% % % colormap pink
% % % axval = axis;
% % % % Add lines for the limiting values of eta2(lambda1)
% % % lambdalimit = linspace(0.5, 1, 51);
% % % etalimit    = 0.5 * (lambdalimit - lambdalimit.^2);
% % % plot3(lambdalimit, zeros(size(etalimit)),  etalimit, 'k-', ...
% % %       lambdalimit, zeros(size(etalimit)), -etalimit, 'k-', ...
% % %       'LineWidth', 1.5)
% % % axis(axval)
% % % axis([0.5 1 0 1.5 -0.002 0.01])
% % % 
