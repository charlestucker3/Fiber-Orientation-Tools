function [t, psi, phi, varargout] = solvePsi2D(L, CI, xi, tEnd, dtIn, ...
                                               nnod, Co, varargin)
%[T, PSI, PHI] = SOLVEPSI2D(L, CI, XI, TEND, DTIN, NNOD, CO) solves for 
%    the planar orientation distribution function PSI(PHI) at times
%    0:DTIN:TEND, for velocity gradient L, interaction coefficient CI, and
%    particle shape factor XI. L can be 2x2 or 3x3.  The initial condition
%    is random-in-plane, unless otherwise specified via VARGIN (see below).
%    The results are returned at NNOD nodal positions spanning 
%    [-pi/2, pi/2).  CO is the maximum Courant number, and may be used to
%    determine a time sub-step.  CO <= 1 is normal.  The method is finite
%    differences and the time integration is fully implicit.
%    T is a vector of the solution times, PHI is a vector of NNOD angles
%    for the nodes, and PSI(I,J) is the distribution function at node I for
%    time T(J).
%
%    If TEND <= 0 then a single steady-state PSI is returned and DTIN is
%    ignored.
%
%[T, PSI, PHI, AV] = SOLVEPSI2D(... also returns the orientation tensor 
%    components.  AV(J,K) is the contracted tensor component K at time
%    T(J).  The contracted components are A(1,1), A(2,2), A(1,2) for K =
%    1:3, respectively.  
%
%[T, PSI, PHI, AV, NSUB] = SOLVEPSI2D(... also returns NSUB, the number of 
%    time sub-steps used for each DTIN step in the output.
%
%[T, PSI, PHI, AV, NSUB, A4V] = SOLVEPSI2D(... also returns A4V, the
%    fourth-order orientation tensor at each time step.  A4V(J,K) is the
%    contracted tensor component of A4 as in tens2vec4.m.  
%
%[T, PSI, PHI] = SOLVEPSI2D(L, ..., CO, SCHEME) controls the finite
%    difference scheme.  Options for SCHEME are 'central' (central
%    differencing, the default), 'upwind' (first-order upwind differencing;
%    not recommended), and 'power' (power-law upwinding; useful for low-CI
%    solutions on relatively coarse grids).  
%
%[T, PSI, PHI] = SOLVEPSI2D(L, ..., CO, SCHEME, PSIZERO) uses PSIZERO as 
%    the initial condition.  PSIZERO must be a vector of length NNOD.

% -- Parse input arguments and set default values
if tEnd > 0  % It's a transient calculation
    calcType = 'transient';
else % It's a steady-state calculation
    calcType = 'steady';
end

if nargin < 8
    scheme = 'central';           % The default scheme
else
    scheme = lower(varargin{1});  % User-selected scheme
end

% --- Set up grid, velocities, time steps, etc.

% Angle values
dphi = pi/nnod;                            % Node spacing, radians
phi = linspace(-pi/2, (pi/2)-dphi, nnod);  % Nodal phi values, radians

% Angular velocity values. phiDot(i) is halfway between phi(i) and phi(i+1)
% (i.e., phiDot(i) is to the East of psi(i)).
D = 0.5*(L+L');
W = 0.5*(L-L');
gammaDot = sqrt(2*trace(D*D));
phiDot = -W(1,2) - xi*D(1,1)*sin(2*(phi+(dphi/2))) ...
                 + xi*D(1,2)*cos(2*(phi+(dphi/2)));
% Rotary diffusivity
Dr = CI * gammaDot;   

% Time step selection, number of time steps and number of sub-steps
if strcmp(calcType, 'transient')
    % Number of time steps for solution reporting
    nsteps = ceil(tEnd/dtIn);
    % Time step between stored solutions
    dtSol = tEnd/nsteps;
    % The time for each stored solution
    t = linspace(0, tEnd, nsteps+1);
    % Courant-number limit on time step
    dtCourant = Co * dphi / max(abs(phiDot));
    % Number of sub-steps per reporting step, if required
    nsub = ceil(dtSol/dtCourant);
    % Time step size for solution sub-step; used to form K
    dt = dtSol/nsub;
else
    % The steady-state strategy is to take ten huge time steps,
    % which are stable for the implicit scheme.  This avoids
    % potential ill-conditioning of the one-step steady-state equations.
    nsteps = 10; % Set up for 10 steps; only the last will be returned
    dt     = 100 * gammaDot;
    nsub   =  1;  % Placeholder for the output
    t      =  0;  % Placeholder for the output
end

% --- Form "stiffness matrix" K, which is tridiagonal plus entries in the 
%     upper right and lower left corners due to the periodic BCs
K = zeros(nnod, nnod);% Scalar magnitude of the rate of deformation
for i = 1:nnod
    % Indices of nodes to the East and West of node i.
    % This helps handle the periodic boundary condition.
    if i == 1
        iW = nnod;
    else
        iW = i-1;
    end
    if i == nnod
        iE = 1;
    else
        iE = i+1;
    end
    
    % Coefficients for this line of the K matrix
    Fw = phiDot(iW)/dphi;   % Advection, west side
    Fe = phiDot(i) /dphi;   % Advection, east side
    Dw = Dr/(dphi)^2;       % Diffusion, west side
    De = Dr/(dphi)^2;       % Diffusion, east side
    switch scheme
        case 'central'
            aW = Dw + Fw/2;
            aE = De - Fe/2;
        case 'upwind'
            aW = Dw + max( Fw, 0);
            aE = De + max(-Fe, 0);
        case 'power'
            Pew = abs(Fw/Dw);  % Grid Peclet number, west edge
            Pee = abs(Fe/De);  % Grid Peclet number, east edge
            aW = Dw*max(0,(1-0.1*Pew)^5) + max( Fw, 0); 
            aE = De*max(0,(1-0.1*Pee)^5) + max(-Fe, 0); 
        otherwise
            error('%s is not a legal scheme', scheme)
    end
    % For all methods,the central node coefficient is
    aP = aW + aE + Fe - Fw;

    % Add the proper entries to K
    K(i,iW) = -dt*aW;
    K(i,i)  = 1 + dt*aP;
    K(i,iE) = -dt*aE;
end 

% --- Main solution
% Storage for the results
psi  = zeros(nnod, nsteps+1);
% Set the initial condition
if nargin < 9
    % Use a random-in-plane IC
    psi(:,1) = 1/(2*pi);
else
    % Use the user-supplied IC
    psi(:,1) = varargin{2};
end
% Compute K-inverse once, since it doesn't change
Kinv = K\eye(nnod);
% Cycle the time steps and any substeps
for i = 1:nsteps
    psiTmp = psi(:,i);  % Most recent stored solution
    for j = 1:nsub
        % Sub-step the solution
        psiTmp = Kinv * psiTmp;
    end
    % Store the result for this big step
    psi(:,i+1) = psiTmp;
end
        
if strcmp(calcType, 'steady')
    % Return only the final time step
    psi = psi(:,end);
end

% --- Provide any optional output arguments
if nargout >= 4
    % -- Compute the orientation tensor values at each output time
    % Storage for the results
    [~, nvals] = size(psi);  % nvals = number of data sets
    Av = zeros(nvals, 3);
    % Matrix of pp values for each phi; used to convert psi to tensor
    pp = [cos(phi).^2; sin(phi).^2; cos(phi).*sin(phi)];
    for i = 1:nvals
        Av(i,:) = pp * psi(:,i) * (2*dphi);
    end
    varargout{1} = Av;
    
    if nargout >= 6
    % -- Compute the fourth-order orientation tensor at each time,
    %    using the contracted form as in tens2vec4.m
    A4v = zeros(nvals, 5);
    pppp = [cos(phi).^4; sin(phi).^4; cos(phi).^2 .* sin(phi).^2; ...
            cos(phi).^3 .* sin(phi); cos(phi) .* sin(phi).^3];
    for i = 1:nvals
        A4v(i,:) = pppp * psi(:,i) * (2*dphi);
    end
    varargout{3} = A4v;
    end
end

if nargout >= 5
    % Also report the number of sub-steps
    varargout{2} = nsub;
end

return