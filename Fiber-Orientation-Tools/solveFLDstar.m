function [tout, N, Ln, Lw, varargout] = solveFLDstar(Nzero, Lstar, t, ...
                                                     S, varargin)
%[TOUT, N, LN, LW] = FLDSTAR(NZERO, LSTAR, T, S) solves the Phelps-Tucker 
%     fiber length model in dimensionless form, for constant unbreakable
%     fiber length Lub.  Usually this means constant matrix viscosity and
%     shear rate.  Nzero is a column vector giving the initial condition
%     for the number-based fiber length distribution and LSTAR is a vector
%     giving the corresponding fiber lengths.  LSTAR(i) = i*Delta_L/Lub.  T
%     is a vector of times at which the solution should be calculated, and
%     S is the dimensionless parameter controlling the distribution of
%     child lengths for each parent fiber.  TOUT is a vector of actual
%     solution times.
%
%[TOUT, N, LN, LW] = FLDSTAR(NZERO, LSTAR, T, S, PMODEL) uses PMODEL for
%     the breakage rate vs. fiber length model.  
%     Choices for PMODEL are:
%          'Phelps'    The Phelps-Tucker model (the default).  This has
%                      zero breakage rate for L < Lub.
%          'Durin'     The model of Durin et al., Composites A, v48, p47,
%                      2013.  This allows some breakage for L < Lub, and
%                      has full breakage for L > Lub.
%
%     The function returns a matrix N, where N(I,J) is the number-based FLD
%     for length LSTAR(I) at time TOUT(J).  N(:,1) equals the initial
%     condiiton.  Also returned are vectors of the number-average and
%     weight-average fiber length, LN and LW, at times
%     TOUT.  
%
%[TOUT, N, LN, LW, W] = FLDSTAR(NZERO, LSTAR, T, S) also returns a 
%     normalized weight-based FLD at each time in the matrix W, which has
%     the same dimensions as N.

%  Problem size and solution storage
nlen = length(Lstar);   % Number of fiber lengths in the discrete FLD
Nzero = reshape(Nzero, nlen, 1);   % Ensure Nzero is a column vector
Lstar = reshape(Lstar, nlen, 1);   % Ensure Lstar is a column vector

% Parse variable input arguments
if nargin >= 5
    Pmodel = varargin{1};  % Should be 'Phelps' or 'Durin'
else
    Pmodel = 'Phelps';      % The default
end

% Transition matrix for the discrete FLD
% - First, find the scaled matrix, without breakage rates.
Rhat = fldRstar(nlen, S);   % Dimensionless matrix w/o breakage rates

% -- The breakage rates depens on the model used
switch lower(Pmodel)
    case 'phelps'
        P = 1 - exp(1 - Lstar.^4);  % Breakage rate, ...
        P(P<0) = 0;                 % ... but never less than zero
    case 'durin'  % Eqns. 25a) and (25b) in Durin (2013).
        P = (1-exp(-Lstar.^4)) / (1-exp(-1));
        P(P>1) = 1;
    otherwise
        error('Illegal value of Pmodel.')
end
P = reshape(P, 1, nlen);    % Make sure P is a row vector

% Multiply each column of Rhat by the corresponding P 
Rhat = repmat(P, nlen, 1) .* Rhat; 
% Rhat is now the complete matrix so that dN/dt = Rhat * N.

% Runge-Kutta time integration for the solution
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Options for ode45
[tout, N] = ode45(@(t,Nstep) Rhat*Nstep, t, Nzero, options);
N = N';  % Make each column of N the FLD for a time step

% Post-process to get Ln, Lw, w
Ln = zeros(size(tout));
Lw = zeros(size(tout));
w  = zeros(size(N));
for k = 1:length(tout)
    Ln(k)  = sum(N(:,k).*Lstar)    / sum(N(:,k));
    Lw(k)  = sum(N(:,k).*Lstar.^2) / sum(N(:,k).*Lstar);
    if nargout >= 4
        % Weight-based distribution
        w(:,k) = N(:,k).*Lstar;
        w(:,k) = w(:,k) / sum(w(:,k));
    end
end

if nargout >= 4  % Return the weight-based distribution fcn. at each tout
    varargout{1} = w;
end

return
