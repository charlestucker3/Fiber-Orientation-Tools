function [t, Av, varargout] = solvePsiARD(L, diffModel, diffParam, xi, ...
                              tEnd, dtIn, ntheta, nphi, Co, varargin)
%[T, AV] = SOLVEPSIARD(L, DIFFMODEL, DIFFPARAM, XI, TEND, DTIN, NTHETA, 
%    NPHI, CO) solves an anisotropic rotary diffusion (ARD) model for the
%    3-D orientation distribution function PSI(THETA, PHI) at times
%    0:DTIN:TEND, for velocity gradient L and particle shape factor XI. L
%    must be 3x3.
%    DIFFMODEL is a character string specifying the ARD model for the
%    rotary diffusion tensor C.  DIFFPARAM gives the model parameters, and
%    its interpretation is dependent on the DIFFMODEL.  Options for
%    DIFFMODEL and DIFFPARAM are:
%      'F'  Folgar-Tucker isotropic model, for testing. DIFFPARAM(1) = CI
%      'I'  Moldex3D iARD model.   DIFFPARAM = [CI, CM]
%      'P'  Moldex3D pARD model.   DIFFPARAM = [CI, omega]
%      'M'  Moldflow's MRD model.  DIFFPARAM = [CI, D2, D3], with D1 = 1  
%      'W'  Wang model (WPT)       DIFFPARAM = [b1, b3]
%    See Favaloro and Tucker, Composites Part A, 126, 105605 (2019) or
%    Section 5.6 of Fundamentals of Fiber Orientation (C.L. Tucker, Hanser,
%    2022) for more information on the various ARD models.
%
%    The initial condition is isotropic, unless otherwise specified via
%    VARARGIN (see below). The solution is calculated on grid of NTHETA by
%    NPHI cells covering a hemisphere. CO is the maximum Courant number,
%    and may be used to determine a time sub-step.  CO <= 1 is normal.  The
%    method is finite differences and the time integration is fully
%    implicit. T is a vector of the solution times, and AV(K.J) contains
%    the Kth component of the second-order orientation tensor in contracted
%    form for time T(J).
%
%    If TEND <= 0 then a single steady-state result is returned and DTIN is
%    ignored.
%
%[T, AV, A4V] = SOLVEPSIARD(... also returns the fourth-order orientation 
%    tensor components at each time step.  A4V(M,J) contains the Mth
%    component of A4 in vector form, arranged as in TENS2VEC4, at time
%    T(J).
%
%[T, AV, A4V, PSI, THETA, PHI, PX, PY, PZ, NSUB] = SOLVEPSIARD(... also 
%    returns the distribution function PSI at each time in T, at cell
%    centroid locations THETA and PHI.  THETA and PHI are NTHETA x NPHI
%    arrays.  PSI(:,J) is the solution at time T(J), and RESHAPE(PSI(:,J),
%    NTHETA, NPHI) will produce an array that matches THETA and PHI. PX, PY
%    and PZ are the Cartesian components of the orientation vector at each
%    cell centroid, and are useful for SURF plots of PSI. NSUB is the
%    number of time sub-steps used for each DTIN step in the output.
%
%[T, AV] = SOLVEPSIARD(L, ..., CO, SCHEME) controls the finite
%    difference scheme.  Options for SCHEME are 'central' (central
%    differencing, the default), 'upwind' (first-order upwind differencing;
%    not recommended), and 'power' (power-law upwinding; useful for low-CI
%    solutions on relatively coarse grids).  
%
%[T, AV] = SOLVEPSIARD(L, ..., CO, SCHEME, PSIZERO) uses PSIZERO as 
%    the initial condition.  PSIZERO must be a vector of length
%    NTHETA*NPHI.

%%
% --- Parse input arguments and set default values
if tEnd > 0  % It's a transient calculation
    calcType = 'transient';
else % It's a steady-state calculation
    calcType = 'steady';
end

if nargin < 10
    scheme = 'central';           % The default scheme
else
    scheme = lower(varargin{1});  % User-selected scheme
end

% --- Set up grid geometry
%     The solution space is the hemisphere 0 <= theta <= pi, 
%     0 <= phi <= pi, since the solution is periodic.

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
%   North:  thetaNS(i,  j),   phiNS(i,  j)
%   South:  thetaNS(i+1,j),   phiNS(i+1,j)
%   West:   thetaEW(i,  j),   phiEW(i,  j)
%   East:   thetaEW(i,  j+1), phiEW(i,  j+1)

% Cell areas.  These vary only with theta, so a vector suffices.
% Here we use the exact result.
Acell = dphi * (cos(thetaNS(1:end-1,1)) - cos(thetaNS(2:end,1)));

% Indices in the grand "stiffness" matrix K for each (i,j) cell
% ans its eight neighbors.  
% ig(i,j) and jg(i,j) are the gridded values of i and j.
[ig, jg] = ndgrid(1:ntheta, 1:nphi);
kP  = ig + (jg-1)*ntheta;   % Each node's "home" index; the k = f(i,j)eqn.
kN  = 1 + mod(ig-2,ntheta) + (jg-1)*ntheta;  % North neighbor
kS  = 1 + mod(ig,  ntheta) + (jg-1)*ntheta;  % South neighbor
kE  = ig + mod(jg,nphi)*ntheta;              % East  neighbor
kE(:,nphi) = ntheta-ig(:,nphi)+1 + mod(jg(:,nphi),  nphi)*ntheta; % East BC
kW  = ig + mod(jg-2,nphi)*ntheta;            % West  neighbor
kW(:,1)    = ntheta-ig(:,1)   +1 + mod(jg(:,1)-2,   nphi)*ntheta; % West BC
kNE = 1 + mod(ig-2,ntheta) + mod(jg,nphi)*ntheta;   % Northeast neighbor
kNE(:,nphi) = 1 + mod(ntheta-ig(:,nphi)+1,ntheta) ...
            + mod(jg(:,nphi),nphi)*ntheta;          % East BC
kSE = 1 + mod(ig,ntheta)   + mod(jg,nphi)*ntheta;   % Southeast neighbor
kSE(:,nphi) = 1 + mod(ntheta-ig(:,nphi)-1,ntheta) ...
            + mod(jg(:,nphi),nphi)*ntheta;          % East BC
kSW = 1 + mod(ig,ntheta)   + mod(jg-2,nphi)*ntheta; % Southwest neighbor
kSW(:,1)    = 1 + mod(ntheta-ig(:,1)   -1,ntheta) ...
            + mod(jg(:,1)-2,nphi)*ntheta;           % West BC
kNW = 1 + mod(ig-2,ntheta) + mod(jg-2,nphi)*ntheta; % Northwest neighbor
kNW(:,1)    = 1 + mod(ntheta-ig(:,1)   +1,ntheta) ...
            + mod(jg(:,1)-2,nphi)*ntheta;           % West BC

% p-vector components at cell centroids.  px(i,j) is the x components of p
% for cell (i,j).
px = cos(phi).*sin(theta);
py = sin(phi).*sin(theta);
pz =           cos(theta);

% The ppA matrix, used to calculate the second-order orientation tensor
% at each time. Av(m,n) will be the mth contracted component of
% the second-order orientation tensor at time step n.  Av is 6 x nsteps+1.
% ppA4 will be used to calculation the fourth-order orientation tensor. 
% A4v(m,n) will be the mth component of the fourth-order tensor at time
% step n.  The order of components matches tens2vec4.
nnod = ntheta *nphi;   % Number of cells (nodes) where psi is calculated
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
% Now the tensor calculations will be Av = 2 * ppA * psi
% and A4v = 2 * ppA4 * psi
% The factors of 2 appear because the solution grid only covers half the
% unit sphere.

% --- Jeffery velocities: vtheta = d(theta)/dt at the NS cell boundaries
%     and vphi = sin(theta)*d(phi)/dt at the EW boundaries.
D = (L+L')/2;
W = (L-L')/2;
gammaDot = sqrt(2*trace(D*D));
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



% ---  Storage for the results. 
% Psi for each time step is a column vector, 
% indexed by k = i + (j-1)*ntheta.
% psiGrid = reshape(psi(:,n), ntheta, nphi) will give psi in an (i,j) grid.
psi  = zeros(nnod, nsteps+1);
Av   = zeros(6,    nsteps+1);  % Store second-order orientation tensors
A4v  = zeros(15,   nsteps+1);  % Store fourth-order orientation tensors

% --- Initial condition
if nargin < 11
    % Use a isotropic IC
    psi(:,1) = 1/(4*pi);
else
    % Use the user-supplied IC
    psi(:,1) = varargin{2};
end
% Orientation tensors at the initial condition
Av( :,1) = 2 * ppA  * psi(:,1);
A4v(:,1) = 2 * ppA4 * psi(:,1);

%% --- Main solution loop
for n = 1:nsteps        % Loop over stored time steps
    psiTmp = psi(:,n);  % Most recent stored solution
    for m = 1:nsub      % Loop over sub-steps
        
        % --- Calculate the orientation tensor A and the diffusivity tensor
        %     C at the start of the sub-step
        A = vec2tens(2 * ppA * psiTmp);   % A, in 3x3 matrix form
        
        switch upper(diffModel)
            case 'F'  % Folgar-Tucker, isotropic rotary diffusion
                C = diffParam(1) * eye(3);
            case 'I'  % Moldex3D iARD model
                C = diffParam(1) * (eye(3) - 4*diffParam(2)*D*D/(gammaDot^2));
            case 'P'  % Moldex3D pARD model
                [~, R] = eigsort(A);  % Sorted eigenvectors of A
                C = diffParam(1) * (R(:,1)*R(:,1)' + diffParam(2)*R(:,2)*R(:,2)' ...
                    + (1-diffParam(2))*R(:,3)*R(:,3)');
            case 'M'  % Moldflow MRD model
                [~, R] = eigsort(A);  % Sorted eigenvectors of A
                C = diffParam(1) * (R(:,1)*R(:,1)' + diffParam(2)*R(:,2)*R(:,2)' ...
                    + diffParam(3)*R(:,3)*R(:,3)');
            case 'W'  % Wang/WPT/2-constant model
                C = diffParam(1)*eye(3) + diffParam(2)*A*A;
            otherwise
                error('%s is not a legal value of model', diffModel)
        end
        
        % -- Theta and phi components of surface diffusivities, at 
        %    the midpoints of the cell boundaries.
        %  - DnNS and DtNS are the normal (theta-theta) and tangential
        %    (theta-phi) components at the north/south cell boundaries.
        DnNS = zeros(size(thetaNS));
        DtNS = zeros(size(thetaNS));
        for i = 1:ntheta+1
            for j = 1:nphi
                % T is the tangent operator at this point.  It transforms a
                % space tensor like C to a surface tensor with (theta,phi)
                % components.
                T = [cos(thetaNS(i,j))*cos(phiNS(i,j)), -sin(phiNS(i,j)); ...
                     cos(thetaNS(i,j))*sin(phiNS(i,j)),  cos(phiNS(i,j)); ...
                    -sin(thetaNS(i,j)),                  0];
                % Dsurf is the local surface diffusivity, 2x2 form, with
                % (theta,phi) components.
                Dsurf = gammaDot * T'*C*T;
                DnNS(i,j) = Dsurf(1,1);  % = Dsurf(theta,theta)
                DtNS(i,j) = Dsurf(1,2);  % = Dsurf(theta,phi)
            end
        end
        % -  DnEW and DtEW are the normal (phi-phi) and tangential
        %    (phi-theta) components at the east/west cell boundaries.
        DnEW = zeros(size(thetaEW));
        DtEW = zeros(size(thetaEW));
        for i = 1:ntheta
            for j = 1:nphi+1
                % T is the tangent operator at this point.  
                T = [cos(thetaEW(i,j))*cos(phiEW(i,j)), -sin(phiEW(i,j)); ...
                     cos(thetaEW(i,j))*sin(phiEW(i,j)),  cos(phiEW(i,j)); ...
                    -sin(thetaEW(i,j)),                  0];
                % Dsurf is the local surface diffusivity, 2x2 form, with
                % (theta,phi) components.
                Dsurf = gammaDot * T'*C*T;
                DnEW(i,j) = Dsurf(2,2);  % = Dsurf(phi,phi)
                DtEW(i,j) = Dsurf(2,1);  % = Dsurf(phi,theta)
            end
        end
        
        %% --- Form the "stiffness" matrix K in sparse form
        %      The finite difference equations have the form
        %      K*psi(:,n) = psi(:,n-1)
        %  The efficient way to form a large sparse K in Matlab is to build vectors
        %  of the indices where the non-zero entries will go, and the values of
        %  those entries, and only then form the sparse-matrix structure.  This
        %  scheme is faster than directly assigning entries to K for ntheta >= 90.
        krow = zeros(9*nnod,1);    % Row indices in K for non-zero entries
        kcol = zeros(size(krow));  % Column indices in K for non-zero entries
        kval = zeros(size(krow));  % Values of entries in K that don't change
        
        %  The loops for forming K are over the (i,j) grid indices
        k = 0;  % k will index the entries in krow, kcol and kval
        for i = 1:ntheta
            for j = 1:nphi
                
                % -- Advection coefficients for this line of the K matrix
                % -  Advection terms, north/south/east/west edges
                Fn = vtheta(i,j)   * sin(thetaNS(i,j))   * dphi;
                Fs = vtheta(i+1,j) * sin(thetaNS(i+1,j)) * dphi;
                Fw = vphi(i,j)   * dtheta;
                Fe = vphi(i,j+1) * dtheta;
                
                % -  Normal-direction diffusion terms, n/s/e/w edges
                Dn =  DnNS(i,j)   * sin(thetaNS(i,j))   * dphi/dtheta;
                Ds =  DnNS(i+1,j) * sin(thetaNS(i+1,j)) * dphi/dtheta;
                Dw = (DnEW(i,j)  /sin(theta(i,j))) * dtheta/dphi;
                De = (DnEW(i,j+1)/sin(theta(i,j))) * dtheta/dphi;
                
                % -  Tangential-direction diffusion terms, n/s/e/w edges
                Dnt = DtNS(i,  j)/4;
                Dst = DtNS(i+1,j)/4;
                Dwt = DtEW(i,  j)/4;
                Det = DtEW(i,j+1)/4;
                
                % -- Node coefficients aN, etc. depend on upwinding scheme
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
                
                % Tangential-diffusion coefficients
                bN  =  Dwt - Det;
                bS  = -bN;
                bW  =  Dnt - Dst;
                bE  = -bW;
                bNE = -Dnt - Det;
                bSE =  Dst + Det;
                bNW =  Dnt + Dwt;
                bSW = -Dst - Dwt;
                
                % Add the proper entries to krow, kcol, kval.
                % Note that kP, kN, kS, etc. were built earlier
                krow(k+(1:9)) = kP(i,j);  % All entries are on row kP
                kcol(k+1) = kP(i,j);      % The cell-center node
                kval(k+1) = 1 +   aP * dt/Acell(i);
                kcol(k+2) = kN(i,j);      % North neighbor
                kval(k+2) = -(aN+bN) * dt/Acell(i);
                kcol(k+3) = kS(i,j);      % South neighbor
                kval(k+3) = -(aS+bS) * dt/Acell(i);
                kcol(k+4) = kE(i,j);      % East neighbor
                kval(k+4) = -(aE+bE) * dt/Acell(i);
                kcol(k+5) = kW(i,j);      % West neighbor
                kval(k+5) = -(aW+bW) * dt/Acell(i);
                kcol(k+6) = kNE(i,j);     % Northeast neighbor
                kval(k+6) = -bNE     * dt/Acell(i);
                kcol(k+7) = kSE(i,j);     % Southeast neighbor
                kval(k+7) = -bSE     * dt/Acell(i);
                kcol(k+8) = kSW(i,j);     % Southwest neighbor
                kval(k+8) = -bSW     * dt/Acell(i);
                kcol(k+9) = kNW(i,j);     % Northwest neighbor
                kval(k+9) = -bNW     * dt/Acell(i); 
                
                k = k+9; % Increment the m counter; there are 9 entries per cell
            end
        end  % krow, kcol and kval are now fully formed
        % Form the sparse matrix K
        K = sparse(krow, kcol, kval, nnod, nnod);
        
        % -- Solve the sub-step
        psiTmp = K \ psiTmp;
    end  % End of the sub-step solution loop
    
    % -- Store the results for this big step
    psi(:,n+1) = psiTmp;
    Av(:, n+1) = 2 * ppA  * psiTmp;
    A4v(:,n+1) = 2 * ppA4 * psiTmp;
    
end % End of the loop over stored time steps


        
if strcmp(calcType, 'steady')
    % Return only the final time step
    psi = psi(:,end);
    Av  = Av( :,end);
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
    % Detailed distribution function results
    varargout{8} = nsub;
end


return