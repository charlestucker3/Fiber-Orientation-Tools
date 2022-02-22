function [t, Av, varargout] = solvePsi3D(L, CI, xi, tEnd, dtIn, ...
                                         ntheta, nphi, Co, varargin)
%[T, AV] = SOLVEPSI3D(L, CI, XI, TEND, DTIN, NTHETA, NPHI, CO) solves for 
%    the 3-D orientation distribution function PSI(THETA, PHI) at times
%    0:DTIN:TEND, for velocity gradient L, interaction coefficient CI, and
%    particle shape factor XI. L must be 3x3.  The initial condition
%    is isotropic, unless otherwise specified via VARARGIN (see below).
%    The solution is calculated on grid of NTHETA by NPHI cells covering a
%    hemisphere. CO is the maximum Courant number, and may be used to
%    determine a time sub-step.  CO <= 1 is normal.  The method is finite
%    differences and the time integration is fully implicit.
%    T is a vector of the solution times, and AV(K,J) contains the Kth
%    component of the second-order orientation tensor in contracted form
%    for time T(J).  
%
%    If TEND <= 0 then a single steady-state result is returned and DTIN is
%    ignored.
%
%[T, AV, A4V, PSI, THETA, PHI] = SOLVEPSI3D(...) also 
%    returns the fourth-order tensor in A4V and the distribution function
%    PSI at each time in T.  A4V(M,J) contains the Mth component of A4 in
%    vector form, arranged as in TENS2VEC4, at time T(J).  PSI(:,J) gives
%    the solution at time T(J) at cell centroid locations THETA and PHI.
%    THETA and PHI are NTHETA x NPHI arrays.  RESHAPE(PSI(:,J), NTHETA,
%    NPHI) produces an array that matches THETA and PHI. 
%
%[T, AV, A4V, PSI, THETA, PHI, PX, PY, PZ, NSUB] = SOLVEPSI3D(...) also 
%    returns, PX, PY and PZ, the Cartesian components of the orientation
%    vector at each cell centroid.  These are useful for SURF plots of PSI.
%    NSUB is the number of time sub-steps used for each DTIN step in the
%    output.  Items may be omitted from the end of the output list.
%
%[T, AV] = SOLVEPSI3D(L, ..., CO, SCHEME) controls the finite
%    difference scheme.  Options for SCHEME are 'central' (central
%    differencing, the default), 'upwind' (first-order upwind differencing;
%    not recommended), and 'power' (power-law upwinding; useful for low-CI
%    solutions on relatively coarse grids).  
%
%[T, AV] = SOLVEPSI3D(L, ..., CO, SCHEME, PSIZERO) uses PSIZERO as 
%    the initial condition.  PSIZERO must be a vector of length
%    NTHETA x NPHI.  


%%
% --- Parse input arguments and set default values
if tEnd > 0  % It's a transient calculation
    calcType = 'transient';
else % It's a steady-state calculation
    calcType = 'steady';
end

if nargin < 9
    scheme = 'central';           % The default scheme
else
    scheme = lower(varargin{1});  % User-selected scheme
end

% --- Set up grid geometry
%     The solution space is the hemisphere 0 <= theta <= pi, 
%     0 <= phi <= pi, since the solution is periodic.
%     The grid directions are defined as on a globe:
%       north is the direction of decreasing theta,
%       south is the direction in increasing theta;
%       east  is the direction of decreasing phi;
%       west  is the direction of increasing phi.

% Angle increments
dtheta = pi/ntheta;
dphi   = pi/nphi;

% Angle coordinates at cell centroids.  
[theta,phi] = ndgrid(linspace(dtheta/2, pi-dtheta/2, ntheta), ...
                     linspace(dphi/2,   pi-dphi/2,   nphi)  );

% Angle coordinates at cell boundaries, where these differ from centroids
%  North/south cell boundaries have different theta values
[thetaNS, phiNS] = ndgrid(linspace(0,      pi,          ntheta+1), ...
                          linspace(dphi/2, pi-dphi/2,   nphi)    );
%  East/west cell boundaries have different theta values
[thetaEW, phiEW] = ndgrid(linspace(dtheta/2, pi-dtheta/2, ntheta), ...
                          linspace(0,        pi,          nphi+1));
% Now, in grid indices, cell (i,j) is located at theta(i,j), phi(i,j),
% and its boundaries are located at:
%   North  thetaNS(i,  j),   phiNS(i,  j)
%   South  thetaNS(i+1,j),   phiNS(i+1,j)
%   West   thetaEW(i,  j),   phiEW(i,  j)
%   East   thetaEW(i,  j+1), phiEW(i,  j+1)

% Cell areas.  These vary only with theta, so a vector suffices.
% Acell(i) will be the area of cell (i,j).  
% Here we use the exact result.
Acell = dphi * (cos(thetaNS(1:end-1,1)) - cos(thetaNS(2:end,1)));

% p-vector components at cell centroids.  px(i,j) contains the x components
% of p for cell (i,j).
px = cos(phi).*sin(theta);
py = sin(phi).*sin(theta);
pz =           cos(theta);

% --- Jeffery velocities: vtheta = d(theta)/dt at the NS cell boundaries
%     and vphi = sin(theta)*d(phi)/dt at the EW boundaries.
D = (L+L')/2;
W = (L-L')/2;
gammaDot = sqrt(2*trace(D*D));
Dr = CI * gammaDot;   
% North/South boundary velocities, done in a loop for clarity
vtheta = zeros(size(thetaNS));
for i = 1:ntheta+1
    for j = 1:nphi
        % The local p, in column vector form
        p = [cos(phiNS(i,j)) * sin(thetaNS(i,j)); ...
             sin(phiNS(i,j)) * sin(thetaNS(i,j)); ...
                             cos(thetaNS(i,j))];
        % dp/dt from Jeffery's equation
        dpdt = W*p + xi*(D*p - (p'*D*p)*p);
        % theta-direction component of dpdt
        vtheta(i,j) = dpdt(1)*cos(thetaNS(i,j))*cos(phiNS(i,j)) ...
                    + dpdt(2)*cos(thetaNS(i,j))*sin(phiNS(i,j)) ...
                    - dpdt(3)*sin(thetaNS(i,j));
    end
end
% East/West boundary velocities, done in a loop for clarity
vphi = zeros(size(phiEW));
for i = 1:ntheta
    for j = 1:nphi+1
        % The local p, in column vector form
        p = [cos(phiEW(i,j)) * sin(thetaEW(i,j)); ...
             sin(phiEW(i,j)) * sin(thetaEW(i,j)); ...
                             cos(thetaEW(i,j))];
        % dp/dt from Jeffery's equation
        dpdt = W*p + xi*(D*p - (p'*D*p)*p);
        % phi-direction component of dp/dt
        vphi(i,j)   = -dpdt(1)*sin(phiEW(i,j)) + dpdt(2)*cos(phiEW(i,j));
    end
end


% --- Time step selection, number of time steps and number of sub-steps
if strcmp(calcType, 'transient')
    % Number of time steps for solution reporting
    nsteps = ceil(tEnd/dtIn);
    % Time step between stored solutions
    dtSol = tEnd/nsteps;
    % The time for each stored solution
    t = linspace(0, tEnd, nsteps+1);
    % Courant-number limit on time step
    dtCoTheta = Co * dtheta / max(max(abs(vtheta)));
    dtCoPhi   = Co * dphi   / max(max(abs(vphi)));
    dtCourant = min(dtCoTheta, dtCoPhi);
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

%% --- Form "stiffness" matrix K as a sparse matrix.  
%      The finite difference equations have the form 
%      K*psi(:,n) = psi(:,n-1)
nnod = ntheta *nphi;   % Number of cells (nodes) where psi is calculated
%  The efficient way to form a large sparse K in Matlab is to build vectors
%  of the indices where the non-zero entries will go, and the values of
%  those entries, and only then form the sparse-matrix structure.  This
%  scheme is faster than directly assigning entries to K for ntheta >= 90.
krow = zeros(5*nnod,1);    % Row indices in K for non-zero entries
kcol = zeros(size(krow));  % Column indices in K for non-zero entries
kval = zeros(size(krow));  % Values of non-zero entries in K

%  The main loop is over the (i,j) grid indices
m = 0;  % m will index the entries in krow, kcol and kval
for i = 1:ntheta
    for j = 1:nphi
        % Grid indices of the cells to the North and South of (i,j).
        % This implements a boundary condition used by Bay, though I
        % think there is no flux on these zero-area boundaries.  This
        % at least avooids special cases in the assembly of K.
        if i == 1
            iN = ntheta;
            jN = j;
        else
            iN = i-1;
            jN = j;
        end
        if i == ntheta
            iS = 1;
            jS = j;
        else
            iS = i+1;
            jS = j;
        end
        % Grid indices of cells to the East and West of cell (i,j).
        % This handles the periodic boundary condition.
        if j == 1
            iW = ntheta-i+1;
            jW = nphi;
        else
            iW = i;
            jW = j-1;
        end
        if j == nphi
            iE = ntheta-i+1;
            jE = 1;
        else
            iE = i;
            jE = j+1;
        end
        
        % Coefficients for this line of the K matrix
        % - Advection terms, north/south/east/west edges
        Fn = vtheta(i,j)   * sin(thetaNS(i,j))   * dphi; 
        Fs = vtheta(i+1,j) * sin(thetaNS(i+1,j)) * dphi;
        Fw = vphi(i,j)   * dtheta;
        Fe = vphi(i,j+1) * dtheta;
        % - Diffusion terms, north/south/east/west edges (Dr = diffusivity)
        Dn = Dr * sin(thetaNS(i,j))   * dphi/dtheta;
        Ds = Dr * sin(thetaNS(i+1,j)) * dphi/dtheta;
        Dw = (Dr/sin(theta(i,j))) * dtheta/dphi;
        De = (Dr/sin(theta(i,j))) * dtheta/dphi;
        switch scheme
            case 'central'
                aN = Dn + Fn/2;
                aS = Ds - Fs/2;
                aW = Dw + Fw/2;
                aE = De - Fe/2;
            case 'upwind'
                aN = Dn + max( Fn, 0);
                aS = Ds + max(-Fs, 0);
                aW = Dw + max( Fw, 0);
                aE = De + max(-Fe, 0);
            case 'power'
                % Grid Peclet numbers, north/south/east/west edges
                Pen = abs(Fn/Dn);
                Pes = abs(Fs/Ds);
                Pew = abs(Fw/Dw);  
                Pee = abs(Fe/De);  
                % Nodal coefficients, north/south/east/west edges
                aN = Dn*max(0,(1-0.1*Pen)^5) + max( Fn, 0);
                aS = Ds*max(0,(1-0.1*Pes)^5) + max(-Fs, 0);
                aW = Dw*max(0,(1-0.1*Pew)^5) + max( Fw, 0);
                aE = De*max(0,(1-0.1*Pee)^5) + max(-Fe, 0);
            otherwise
                error('%s is not a legal scheme', scheme)
        end
        % For all methods,the central node coefficient is
        aP = aN + aS + aW + aE - Fn + Fs - Fw + Fe;
        
        % Compute the indices in K for each (i,j) used by this cell.
        kP = i  + (j-1) *ntheta;  % Cell center
        kN = iN + (jN-1)*ntheta;  % North edge
        kS = iS + (jS-1)*ntheta;  % South edge
        kE = iE + (jE-1)*ntheta;  % East  edge
        kW = iW + (jW-1)*ntheta;  % West  edge
        % Add the proper entries to krow, kcol, kval
        krow(m+(1:5)) = kP;  % All entries for this cell are on row kP
        kcol(m+1) = kP;      % The cell center
        kval(m+1) = 1 + aP*dt/Acell(i); 
        kcol(m+2) = kN;      % North edge
        kval(m+2) =   - aN*dt/Acell(i);
        kcol(m+3) = kS;      % South edge
        kval(m+3) =   - aS*dt/Acell(i);
        kcol(m+4) = kE;      % East edge
        kval(m+4) =   - aE*dt/Acell(i);
        kcol(m+5) = kW;      % West edge
        kval(m+5) =   - aW*dt/Acell(i);
        
        m = m+5; % Increment the m counter; there are 5 entries per cell 
    end
end  % krow, kcol and kval are now fully formed
% Form the sparse matrix K
K = sparse(krow, kcol, kval, nnod, nnod);

%% --- Main solution
% Storage for the results. 
% Psi for each time step is a column vector, 
% indexed by k = i + (j-1)*ntheta.
% psiGrid = reshape(psi(:,n), ntheta, nphi) will give psi in an (i,j) grid.
psi  = zeros(nnod, nsteps+1);

% Set the initial condition
if nargin < 10
    % Use a isotropic IC
    psi(:,1) = 1/(4*pi);
else
    % Use the user-supplied IC
    psi(:,1) = varargin{2};
end

% Do an LU decomposition of K, since we'll solve many right-hand sides
Klu = decomposition(K, 'lu');
% Cycle the time steps and any substeps
for n = 1:nsteps
    psiTmp = psi(:,n);  % Most recent stored solution
    for m = 1:nsub
        % Sub-step the solution
        psiTmp = Klu \ psiTmp;
    end
    % Store the result for this big step
    psi(:,n+1) = psiTmp;
end

%% Calculate the orientation tensor components at each time.
% Av(m,n) will be the mth contracted component of the second-order orientation
% tensor at time step n.  Av is 6 x nsteps+1.  
% A4v(m,n) will be the mth component of the fourth-order tensor at time
% step n.  The order of components matches tens2vec4.
% Set up a matrix ppA that has products of p components and cell areas,
% such that multiplying ppA * psi will give Av.
% ppA4 is the corresponding matrix for A4v.
px = reshape(px, nnod, 1);  % Arrange p components as column vectors
py = reshape(py, nnod, 1);
pz = reshape(pz, nnod, 1);
Aall = repmat(Acell, nphi, 1); % Aall has cell areas for every cell
ppA = [px.*px.*Aall, py.*py.*Aall, pz.*pz.*Aall, ...
       py.*pz.*Aall, pz.*px.*Aall, px.*py.*Aall]';  % Note the transpose
ppA4 = [px.*px.*px.*px.*Aall, ...
        py.*py.*py.*py.*Aall, ...
        pz.*pz.*pz.*pz.*Aall, ...
        py.*py.*pz.*pz.*Aall, ...
        pz.*pz.*px.*px.*Aall, ...
        px.*px.*py.*py.*Aall, ...
        px.*px.*py.*pz.*Aall, ...
        px.*px.*pz.*px.*Aall, ...
        px.*px.*px.*py.*Aall, ...
        py.*py.*py.*pz.*Aall, ...
        py.*py.*pz.*px.*Aall, ...
        py.*py.*px.*py.*Aall, ...
        pz.*pz.*py.*pz.*Aall, ...
        pz.*pz.*pz.*px.*Aall, ...
        pz.*pz.*px.*py.*Aall]';
% Now the tensor calculation is easy.  The factor of 2 appears because 
% the solution grid only covers half the unit sphere.
Av  = 2 * ppA  * psi;
A4v = 2 * ppA4 * psi;
        
if strcmp(calcType, 'steady')
    % Return only the final time step
    psi = psi(:,end);
    Av  = Av(:,end);
    A4v = A4v(:,end);
end

%% --- Provide any optional output arguments
if nargout >= 3
    % Fourth-order orientation tensor
    varargout{1} = A4v;
end

if nargout >= 4
    % Detailed distribution function results
    varargout{2} = psi;
end

if nargout >= 5
    % Theta values, in grid form
    varargout{3} = theta;
end

if nargout >= 6
    % Phi values in grid form
    varargout{4} = phi;
end

if nargout >= 7
    % px values; change back to grid form
    varargout{5} = reshape(px, ntheta, nphi);
end

if nargout >= 8
    % py values; change back to grid form
    varargout{6} = reshape(py, ntheta, nphi);
end

if nargout >= 9
    % pz values; change back to grid form
    varargout{7} = reshape(pz, ntheta, nphi);
end

if nargout >= 10
    % Number of sub-steps for each main time step
    varargout{8} = nsub;
end


return