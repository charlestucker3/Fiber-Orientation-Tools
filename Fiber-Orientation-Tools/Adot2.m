function Aderiv = Adot2(Av, L, xi, closure, diffModel, varargin)
%ADERIV = ADOT2(AV, L, XI, CLOSURE, DIFFMODEL, DIFFPARAM) returns the 
%     time derivative ADERIV of the fiber orientation tensor for a given
%     velocity gradient tensor L (3x3) and particle shape factor XI, using
%     the closure approximation specified by CLOSURE. DIFFMODEL and
%     DIFFPARAM specify the model and parameters for the rotary diffusion
%     model (see below).  
%
%ADERIV = ADOT2(AV, L, XI, CLOSURE, ...
%              DIFFMODEL, DIFFPARAM, KINMODEL, KINPARAM)
%     Also specifies and orientation kinetics model KINMODEL and associated
%     parameters.  The default is standard (Jeffery) kinetics.  
%
%     AV and ADERIV are 6x1 column vectors, containing the tensor
%     components in contracted notation.  
%
%     CLOSURE is single character that is passed to the closeA4.m.  For
%     example, 'Q' gives the quadratic closure, 'N' the orthotropic version
%     of the natural closure (ORN or ORE), and 'I' the Chung & Kwon IBOF
%     closure.
%
%     DIFFMODEL can be any of the following characters:
%       'J'   Jeffery model; no rotary diffusion.  If DIFFPARAM is present
%             it is ignored.
%       'F'   Folgar-Tucker model.  DIFFPARAM(1) = the interaction
%             coefficient CI.
%       'T'   Phelps-Tucker five-parameter ARD model.  DIFFPARAM(1:5)
%             contains b1:b5.  
%       'I'   iARD model of Tseng et al., 2016.  DIFFPARAM(1) = CI and
%             DIFFPARAM(2) = CM.
%       'P'   pARD model of Tseng et al., 2019. 
%             DIFFPARAM(1) = CI and DIFFPARAM(2) = Omega.
%       'M'   Moldflow Rotational Diffusion (MRD) model.  
%             DIFFPARAM(1) = CI, DIFFPARAM(2) = D2, DIFFPARAM(3) = D3.
%             The code assumes D1 = 1.
%       'W'   Wang, or WPT, model.  DIFFPARAM(1) = b1, DIFFPARAM(2) = b2.  
%
%     KINMODEL can be any of the following:
%       'J'   Standard or Jeffery kinetics.  This is the default if
%             KINMODEL is not specified.  KINPARAM is ignored.
%       'R'   RSC (Reduced Strain Closure) model of Wang, O'Gara and
%             Tucker.  KINPARAM(1) = kappa.
%       'P'   RPR (Retarding Principal Rate) model of Tseng et al., 2013.
%             KINPARAM(1) = kappa = (1 - alpha). 
%       'S'   Slip or SRF model.  The entire expression for dA/dt is
%             multiplied by the scalar factor KINPARAM(1) = kappa.  This
%             model is *NOT* objective, but is included here for comparison
%             since it has appeared frequently in the literature.  

% Unpack variable input arguments
if nargin >= 6
    diffParam = varargin{1};
end
if nargin >= 7
    kinModel = varargin{2};
end
if nargin >= 8
    kinParam = varargin{3};
end

% Preliminary quantities common to all models
A  = vec2tens(Av);                % A in 3x3 matrix form
A4 = closeA4(A, closure);         % 4th-order orientation tensor (6x6)
R4 = diag( [1 1 1 2 2 2], 0);     % used for contracted-notation products
Id4 = diag( [1 1 1 0.5 0.5 0.5], 0); % contracted 4th-order identity tensor
W  = (L-L')/2;                    % vorticity tensor
D  = (L+L')/2;                    % deformation rate tensor
gammadot = sqrt(2*trace(D*D));    % scalar strain rate

% Find sorted eigenvalues and eigenvectors of A.  Not used for every model,
% so this could be made conditional on the model for efficiency.
[Evals, R] = eigsort(A);
                                         
% Quantities dependent on the rotary diffusion model: effective values of
% CI, D, and the target orientation tensor Ahat for rotary diffusion.
switch upper(diffModel)
    case {'J'}   % Jeffery model; no rotary diffusion
        CIhat = 0;
        Ahat  = eye(3)/3;  % Target diffusion tensor; placeholder
        Dhat  = D;
        
    case {'F'}                 % Folgar-Tucker isotropic rotary diffusion
        CIhat = diffParam(1);  % Interaction coefficient
        Ahat  = eye(3)/3;      % Target diffusion tensor is isotropic
        Dhat  = D;
        
    otherwise    % Using an anisotropic rotary diffusion (ARD) model
        % Find the rotary diffusion tensor C for the given model
        switch upper(diffModel)
            
            case 'T'  % Phelps-Tucker 5-constant model
                C = diffParam(1) * eye(3) + diffParam(2)*A ...
                  + diffParam(3)*A*A + diffParam(4)*D*D/gammadot ...
                  + diffParam(5)*D*D/(gammadot^2);
              
            case 'I'  % iARD model
                C = diffParam(1) * ...
                    (eye(3) - 4*diffParam(2)*D*D/(gammadot^2));
                
            case 'P'  % pARD model 
                C = diffParam(1) * ...
                    R * diag([1, diffParam(2), 1-diffParam(2)]) * R';
                
            case 'M'  % MRD model
                C = diffParam(1) * ...
                    R * diag([1, diffParam(2), diffParam(3)]) * R';
                
            case 'W'  % WPT model
                C = diffParam(1)*eye(3) + diffParam(2) *A*A;
                
        end
        % For any ARD model, find the effective values of CI and D
        % and the target orientation tensor Ahat.  
        CIhat = trace(C)/3;   
        Ahat = C/trace(C); 
        Dhat = xi*D - 15*CIhat*gammadot*Ahat;
        xi = 1;  % Use this in dA/dt eqn. below, since it is already 
                 % incorporated into Dhat.  
end  % End of section for handling rotary diffusion models

% Treat orientation kinetics models
if nargin <= 6 % No kinetics model is specified, so use the Jeffery model
    kinModel = 'J'; 
end
switch upper(kinModel)  % Set kappa for use in dAdt eqn., and
                        % adjust A4 and Ahat as needed.
    
    case 'R'  % Reduced Strain Closure (RSC) model
        kappa = kinParam(1);
        % Calculate the fourth-order tensors L4 and M4, 
        % which are based on the eigenvalues and eigenvectors of A,
        % in laboratory coordinates.
        % Start in the principal axes of A.
        L4 = diag([Evals(1) Evals(2) Evals(3) 0 0 0], 0);
        M4 = diag([1 1 1 0 0 0], 0);
        % Then rotate back to laboratory axes
        L4 = rotate4(L4, R);   
        M4 = rotate4(M4, R);   
        
        % We can now modify A4 and Ahat for the RSC model
        A4 = A4 + (1-kappa) * (L4 - M4*R4*A4); 
        Ahat = Ahat - (1-kappa) * vec2tens(M4*R4*tens2vec(Ahat));
        
    otherwise
        kappa = 1;  % Used for all other models in dAdt eqn.
        % For all other models, A4 and Ahat remain unchanged.
        % The slip/SRF and RPR models are implemented after doing a
        % standard-kinetics computation of dA/dt.  
        
end % End of initial modifications for orientation kinetics model.
    

% Main equation for dA/dt, in 3x3 matrix form.  This equation plus the
% calculations above accounts for ARD and RSC, if either or both of those
% models have been chosen.  
% Note that kappa multiplies A in the diffusion term.  
dAdt =  W*A - A*W ...
     + xi * (A*Dhat + Dhat*A - 2*vec2tens(A4*R4*tens2vec(Dhat))) ...
     + 6*CIhat*gammadot * (Ahat - kappa*A); 
 
% Additional modifications for RPR or slip/SRF models, if implemented
switch upper(kinModel)
    
    case 'S'  % Slip/SRF reduces all components of dAdt
        dAdt = kinParam(1)*dAdt;
    
    case 'P'  % RPR model 
        dAdtPrin = R' * dAdt * R;   % dAdt in the principal axes of A
        % Modify the diagonal components of dAdt in these axes
        dAdtRPR  = dAdtPrin - (1-kinParam(1)) * diag(diag(dAdtPrin));
        % Return dAdt to the laboratory axes
        dAdt = R * dAdtRPR * R';
              
    otherwise
        % No action needed   
end
		   
% Return the derivative of A in 6x1 column vector form
Aderiv = tens2vec(dAdt);  
return

