function A4 = closeA4(A2, version)
%A4 = CLOSEA4(A2, VERSION) gives the fourth-order orientation tensor A4
%   in 6x6 contracted notation as a function of the second-order tensor
%   A2, which is in 3x3 form.  VERSION is a single character that 
%   determines the particular closure approximation to use.
%
%   Available values of VERSION are:
%      N  Orthotropic equivalent of the natural closure (ORN)
%      E  Same as N (ORE)
%      I  Fitted IBOF closure of Chung and Kwon (2002)
%      H  Hybrid closure
%      Q  Quadratic closure
%      L  Linear closure (not recommended for flow problems)
%      P  Peaked orthotropic closure (ORP)
%      S  Smooth orthotropic closure (ORS)
%      F  Cintra and Tucker fitted closure (ORF)
%      G  Cintra and Tucker low-CI fit (ORL)
%
%  Note the RSC-type versions of A4 are not available here.  That 
%  must be handled in the function that calculates dA2/dt.

% The function is a giant "switch" statement, based on the version
switch upper(version(1))
    
    case {'Q','L','H'} % --- Quadratic, Linear and Hybrid closures ---
        % Quadratic closure
        Av  = tens2vec(A2);    % makes 6x1 a column vector
        A4Q = Av*Av';          % creates the 6x6 matrix

        % Linear closure
        Iexp = [1 2 3 2 3 1];      % these two arrays expand M to (I,J)
        Jexp = [1 2 3 3 1 2];      % or N to (K,L)
        I2   = eye(3);             % 2nd order identity tensor
        A4L  = zeros(6);
        for m = 1:6
            for n = 1:6
                i = Iexp(m);
                j = Jexp(m);
                k = Iexp(n);
                l = Jexp(n);  % A4L(m,n) will contain the tensor component ijkl
                A4L(m,n) = -(1/35)* ...
                    (I2(i,j)*I2(k,l) + I2(i,k)*I2(j,l) + I2(i,l)*I2(j,k) ) ...
                    + (1/7)* ...
                    ( A2(i,j)*I2(k,l) + A2(i,k)*I2(j,l) + A2(i,l)*I2(j,k) ...
                    + I2(i,j)*A2(k,l) + I2(i,k)*A2(j,l) + I2(i,l)*A2(j,k) );
            end
        end

        % Combine the two forms
        switch upper(version(1))
            case 'Q', f = 1;               % quadratic closure
            case 'L', f = 0;               % linear closure
            case 'H', f = 1 - 27*det(A2);  % hybrid
            otherwise, error('Should not reach this statememt')
        end
        A4 = f*A4Q + (1-f)*A4L;

    case 'I'  % ---- Chung and Kwon IBOF closure ----
        %    Adapted by C. L. Tucker from code provided by Endo Hiroki
        %    of Nagaoka University of Technology, 6/09/20.
        
        % -- Invariants
        AA = A2*A2;  % Will be needed several places
        I2 = 0.5*((trace(A2)^2 - trace(AA)));
        I3 = det(A2);
        
        % -- Fitting parameters
        % -- Pasted in from Chung and Kwon 2002, by CLT.
        a = [0.24940908165786E+02 -0.497217790110754E+00  0.234146291570999E+02;
            -0.435101153160329E+03  0.234980797511405E+02 -0.412048043372534E+03;
            0.372389335663877E+04 -0.391044251397838E+03  0.319553200392089E+04;
            0.703443657916476E+04  0.153965820593506E+03  0.573259594331015E+04;
            0.823995187366106E+06  0.152772950743819E+06 -0.485212803064813E+05;
            -0.133931929894245E+06 -0.213755248785646E+04 -0.605006113515592E+05;
            0.880683515327916E+06 -0.400138947092812E+04 -0.477173740017567E+05;
            -0.991630690741981E+07 -0.185949305922308E+07  0.599066486689836E+07;
            -0.159392396237307E+05  0.296004865275814E+04 -0.110656935176569E+05;
            0.800970026849796E+07  0.247717810054366E+07 -0.460543580680696E+08;
            -0.237010458689252E+07  0.101013983339062E+06  0.203042960322874E+07;
            0.379010599355267E+08  0.732341494213578E+07 -0.556606156734835E+08;
            -0.337010820273821E+08 -0.147919027644202E+08  0.567424911007837E+09;
            0.322219416256417E+05 -0.104092072189767E+05  0.128967058686204E+05;
            -0.257258805870567E+09 -0.635149929624336E+08 -0.152752854956514E+10;
            0.214419090344474E+07 -0.247435106210237E+06 -0.499321746092534E+07;
            -0.449275591851490E+08 -0.902980378929272E+07  0.132124828143333E+09;
            -0.213133920223355E+08  0.724969796807399E+07 -0.162359994620983E+10;
            0.157076702372204E+10  0.487093452892595E+09  0.792526849882218E+10;
            -0.232153488525298E+05  0.138088690964946E+05  0.466767581292985E+04;
            -0.395769398304473E+10 -0.160162178614234E+10 -0.128050778279459E+11]';
        
        
        % -- Beta coefficients
        Beta3 = a(1,1) + a(1,2)*I2 + a(1,3)*I2^2 + a(1,4)*I3 + a(1,5)*I3^2 + ...
            a(1,6)*I2*I3 + a(1,7)*(I2^2)*I3 + a(1,8)*I2*I3^2 + a(1,9)*I2^3 + ...
            a(1,10)*I3^3 + a(1,11)*(I2^3)*I3 + a(1,12)*(I2^2)*(I3^2) + ...
            a(1,13)*I2*I3^3 + a(1,14)*I2^4 + a(1,15)*I3^4 + a(1,16)*(I2^4)*I3 + ...
            a(1,17)*(I2^3)*I3^2 + a(1,18)*(I2^2)*I3^3 + a(1,19)*I2*I3^4 + ...
            a(1,20)*I2^5 + a(1,21)*I3^5;
        
        Beta4 = a(2,1) + a(2,2)*I2 + a(2,3)*I2^2 + a(2,4)*I3 + a(2,5)*I3^2 + ...
            a(2,6)*I2*I3 + a(2,7)*(I2^2)*I3 + a(2,8)*I2*I3^2 + a(2,9)*I2^3 + ...
            a(2,10)*I3^3 + a(2,11)*(I2^3)*I3 + a(2,12)*(I2^2)*(I3^2) + ...
            a(2,13)*I2*I3^3 + a(2,14)*I2^4 + a(2,15)*I3^4 + a(2,16)*(I2^4)*I3 + ...
            a(2,17)*(I2^3)*I3^2 + a(2,18)*(I2^2)*I3^3 + a(2,19)*I2*I3^4 + ...
            a(2,20)*I2^5 + a(2,21)*I3^5;
        
        Beta6 = a(3,1) + a(3,2)*I2 + a(3,3)*I2^2 + a(3,4)*I3 + a(3,5)*I3^2 + ...
            a(3,6)*I2*I3 + a(3,7)*(I2^2)*I3 + a(3,8)*I2*I3^2 + a(3,9)*I2^3 + ...
            a(3,10)*I3^3 + a(3,11)*(I2^3)*I3 + a(3,12)*(I2^2)*(I3^2) + ...
            a(3,13)*I2*I3^3 + a(3,14)*I2^4 + a(3,15)*I3^4 + a(3,16)*(I2^4)*I3 + ...
            a(3,17)*(I2^3)*I3^2 + a(3,18)*(I2^2)*I3^3 + a(3,19)*I2*I3^4 + ...
            a(3,20)*I2^5 + a(3,21)*I3^5;
        
        Beta1 = (3/5)*((-1/7) + (1/5)*Beta3*((1/7) + (4/7)*I2 + (8/3)*I3) + ...
            (-Beta4)*((1/5) - (8/15)*I2 - (14/15)*I3) + ...
            (-Beta6)*((1/35) - (24/105)*I3 - (4/35)*I2 + (16/15)*I2*I3 + ...
            (8/35)*I2^2));
        
        Beta2 = (6/7)*(1 - (1/5)*Beta3*(1+4*I2) + (7/5)*Beta4*((1/6)-I2) - ...
            Beta6*(-(1/5) + (2/3)*I3 + (4/5)*I2 - (8/5)*I2^2));
        
        Beta5 = -(4/5)*Beta3 - (7/5)*Beta4 - (6/5)*Beta6*(1 - (4/3)*I2);
        
        % -- Fully symmetric basis tensors
        Iexp = [1 2 3 2 3 1];  % Arrays to expand the contracted indices
        Jexp = [1 2 3 3 1 2];
        I = eye(3);            % Identity tensor
        S_II   = zeros(6);
        S_IA   = zeros(6);
        S_AA   = zeros(6);
        S_IAA  = zeros(6);
        S_AAA  = zeros(6);
        S_AAAA = zeros(6);
        for m =1:6
            for n =1:6
                i = Iexp(m);
                j = Jexp(m);
                k = Iexp(n);
                l = Jexp(n);
                S_II(m,n) = (1/3)* ...
                    (I(i,j)*I(k,l) + I(i,k)*I(j,l) + I(i,l)*I(j,k) );
                S_IA(m,n) = (1/6)* ...
                    (I(i,j)*A2(k,l) + I(k,l)*A2(i,j) + I(i,k)*A2(j,l) + ...
                    I(j,l)*A2(i,k) + I(i,l)*A2(j,k) + I(j,k)*A2(i,l) );
                S_AA(m,n) = (1/3)* ...
                    (A2(i,j)*A2(k,l) + A2(i,k)*A2(j,l) + A2(i,l)*A2(j,k) );
                S_IAA(m,n) = (1/6)* ...
                    (I(i,j)*AA(k,l) + I(k,l)*AA(i,j) + I(i,k)*AA(j,l) + ...
                    I(j,l)*AA(i,k) + I(i,l)*AA(j,k) + I(j,k)*AA(i,l) );
                S_AAA(m,n) = (1/6)* ...
                    (A2(i,j)*AA(k,l) + A2(k,l)*AA(i,j) + A2(i,k)*AA(j,l) + ...
                    A2(j,l)*AA(i,k) + A2(i,l)*AA(j,k) + A2(j,k)*AA(i,l) );
                S_AAAA(m,n) = (1/6)* ...
                    (AA(i,j)*AA(k,l) + AA(k,l)*AA(i,j) + AA(i,k)*AA(j,l) + ...
                    AA(j,l)*AA(i,k) + AA(i,l)*AA(j,k) + AA(j,k)*AA(i,l) );
            end
        end
                
        % --- Assemble the pieces of the IBOF closure ---
        A4 = Beta1*S_II  + Beta2*S_IA  + Beta3*S_AA + ...
            Beta4*S_IAA + Beta5*S_AAA + Beta6*S_AAAA;

    otherwise   % ---- Form one of the EBOF (orthotropic) closures ----

        %---- Principal values and axes of A2 ----
        [Evals, R] = eigsort(A2);  
        
        %---- Build A4 in principal axes -------
        A4 = zeros(6);
        
        %-- Fill in components on the main diagonal
        switch upper(version(1))
            
            case 'P'   % --- "Peaked" version
                C = [ 0  1  0;
                    0  0  1;
                    1 -1 -1;
                    0  0  0;
                    0  0  0;
                    0  0  0];
                A4 = diag(C * [1, Evals(1), Evals(2)]', 0);
                
            case 'S'   % --- "Smooth" version
                C = [-0.15  1.15 -0.10;
                    -0.15  0.15  0.90;
                    0.60 -0.60 -0.60;
                    0.20 -0.20 -0.20;
                    0.20 -0.20 -0.20;
                    -0.05  0.05  0.30];
                A4 = diag(C * [1, Evals(1), Evals(2)]', 0);
                
            case 'F'    % --- Cintra and Tucker "fitted" version
                C3 = [.060964,  .124711, 1.228982, -.146365, -.082617,  .021654;
                    .371243, -.389402,-2.054116,  .407381,  .646735, -.017978;
                    .555301,  .258844,  .821548, -.262545, -.559003,  .003702;
                    -.369160,  .086169,-2.260574,  .902623,  .357952,  .011208;
                    .318266,  .796080, 1.053907, -.765861, -.288046, -.030219;
                    .371218,  .544992, 1.819756, -.996765, -.822991,  .451773]';
                A4 = diag(C3 * ...
                    [1, Evals(1), Evals(1)^2, Evals(2), Evals(2)^2, Evals(1)*Evals(2)]', 0);
                
            case  'G'   % --- Cintra and Tucker "low-CI" version
                C5 = [.104753,  .162210, 1.288896, -.173176, -.115720,  .010966;
                    .346874, -.451257,-2.187810,  .492970,  .694840, -.041714;
                    .544970,  .286639,  .899635, -.320652, -.578983,  .034013;
                    -.563168, -.028702,-2.402857,  .934196,  .468661,  .094506;
                    .471144,  .864008, 1.133547, -.763206, -.370342, -.100803;
                    .491202,  .652712, 1.975826,-1.068668, -.907158,  .415956]';
                A4 = diag(C5 * ...
                    [1, Evals(1), Evals(1)^2, Evals(2), Evals(2)^2, Evals(1)*Evals(2)]', 0);
                
            case {'E','N'}  % --- equivalent of the Natural closure
                CE = [0.636256796880687,   0.636256796880687,   2.74053289560253;
                    -1.8726629637381  ,  -3.31527229742146 ,  -9.12196509782692;
                    -4.47970873193738 ,  -3.03709939825406 , -12.2570587036254 ;
                    11.9589562332320  ,  11.8273285968852  ,  34.3199018916987 ;
                    3.84459692420086 ,   6.88153952058044 ,  13.8294699121940 ;
                    11.3420924278159  ,   8.43677746778325 ,  25.8684755253884 ;
                    -10.9582626069691  , -15.9120667157641  , -37.7029118029384 ;
                    -20.7277994684132  , -15.1515872606307  , -50.2756431927485 ;
                    -2.11623214471004 ,  -6.48728933641926 , -10.8801761133174 ;
                    -12.3875632855619  ,  -8.63891419284016 , -26.9636915239716 ;
                    9.81598389716748 ,   9.32520343452661 ,  27.3346798054488 ;
                    3.47901510567439 ,   7.74683751713295 ,  15.2650686148651 ;
                    11.7492911177026  ,   7.48146870624441 ,  26.1134914005375 ;
                    0.508041387366637,   2.28476531637958 ,   3.4321384033477 ;
                    4.88366597771489 ,   3.59772251134254 ,  10.6117418066060 ]';
                
                for M = 1:3
                    A4(M,M) = CE(M,1) + CE(M,2)*Evals(1) + CE(M,3)*Evals(2) ...
                        + CE(M,4)*Evals(1)*Evals(2) ...
                        + CE(M,5)*Evals(1)^2  ...
                        + CE(M,6)*Evals(2)^2  ...
                        + CE(M,7)*(Evals(1)^2)*Evals(2) ...
                        + CE(M,8)*Evals(1)*(Evals(2)^2) ...
                        + CE(M,9)*Evals(1)^3  ...
                        + CE(M,10)*Evals(2)^3 ...
                        + CE(M,11)*(Evals(1)^2)*(Evals(2)^2) ...
                        + CE(M,12)*(Evals(1)^3)*Evals(2) ...
                        + CE(M,13)*Evals(1)*(Evals(2)^3) ...
                        + CE(M,14)*Evals(1)^4  ...
                        + CE(M,15)*Evals(2)^4;
                end
                %        -- apply normalization condtitions
                A4(4,4) = 0.5*( 1.0 - 2.0*Evals(1) + A4(1,1) ...
                    - A4(2,2) - A4(3,3) );
                A4(5,5) = 0.5*( 1.0 - 2.0*Evals(2) - A4(1,1) ...
                    + A4(2,2) - A4(3,3) );
                A4(6,6) = 0.5*( -1.0 + 2.0*Evals(1) + 2.0*Evals(2) ...
                    - A4(1,1) - A4(2,2) + A4(3,3) );
                
            case {'B'}  % --- Bingham closure, given by Chaubal and Leal, 
                % Journal of Rheology, v42, pp177-201, 1998.  The rows of
                % the CB coefficient array are the rows of their Table II,
                % p185.
                CB = [0.412251,  0.896044,  -2.02654,   1.71079,   -3.44613,  ...
                    6.13032,  -3.7558,     3.24378,   0.381274,  -3.22331;  ...
                    0.150497, -0.29378,    0.092542,  0.044353,  -0.28173,  ...
                    1.92763,  -1.46408,    0.473044,  0.037776,   0.406202; ...
                    -0.610534,  5.01094,   -8.16518,   3.76426,    4.39396,  ...
                    -6.77247,   2.96854,  -15.1093,   10.6866,     9.92054];
                
                for M = 1:3
                    A4(M,M) =  CB(M,1)  ...
                        + CB(M,2)*Evals(1) ...
                        + CB(M,3)*Evals(1)^2 ...
                        + CB(M,4)*Evals(1)^3 ...
                        + CB(M,5)*Evals(2)  ...
                        + CB(M,6)*Evals(2)^2  ...
                        + CB(M,7)*Evals(2)^3 ...
                        + CB(M,8)*Evals(1)*Evals(2) ...
                        + CB(M,9)*(Evals(1)^2)*Evals(2)  ...
                        + CB(M,10)*Evals(1)*(Evals(2)^2);
                end
                %        -- apply normalization condtitions
                A4(4,4) = 0.5*( 1.0 - 2.0*Evals(1) + A4(1,1) ...
                    - A4(2,2) - A4(3,3) );
                A4(5,5) = 0.5*( 1.0 - 2.0*Evals(2) - A4(1,1) ...
                    + A4(2,2) - A4(3,3) );
                A4(6,6) = 0.5*( -1.0 + 2.0*Evals(1) + 2.0*Evals(2) ...
                    - A4(1,1) - A4(2,2) + A4(3,3) );
                
            otherwise
                fprintf(1,'ortho.m: %c is not a legal value of VERSION\n', version);
                
        end
        
        %-- then fill in remaining components by symmetry
        A4(1,2) = A4(6,6);
        A4(2,1) = A4(1,2);
        A4(2,3) = A4(4,4);
        A4(3,2) = A4(2,3);
        A4(3,1) = A4(5,5);
        A4(1,3) = A4(3,1);
       
        Id4 = diag( [1 1 1 0.5 0.5 0.5], 0);  % contracted 4th-order identity tensor
        R4  = diag( [1 1 1 2   2   2  ], 0);  % inverse of 4th-order identity tensor

        %--- Rotate back to laboratory axes --------
        %---- 6x6 rotation matrices and associated quantities -----
        Qa=[R(1,1)*R(1,1), R(1,2)*R(1,2), R(1,3)*R(1,3), ...
            R(1,2)*R(1,3), R(1,3)*R(1,1), R(1,1)*R(1,2);
            R(2,1)*R(2,1), R(2,2)*R(2,2), R(2,3)*R(2,3), ...
            R(2,2)*R(2,3), R(2,3)*R(2,1), R(2,1)*R(2,2);
            R(3,1)*R(3,1), R(3,2)*R(3,2), R(3,3)*R(3,3), ...
            R(3,2)*R(3,3), R(3,3)*R(3,1), R(3,1)*R(3,2);
            R(2,1)*R(3,1), R(2,2)*R(3,2), R(2,3)*R(3,3), ...
            R(2,2)*R(3,3), R(2,3)*R(3,1), R(2,1)*R(3,2);
            R(3,1)*R(1,1), R(3,2)*R(1,2), R(3,3)*R(1,3), ...
            R(3,2)*R(1,3), R(3,3)*R(1,1), R(3,1)*R(1,2);
            R(1,1)*R(2,1), R(1,2)*R(2,2), R(1,3)*R(2,3), ...
            R(1,2)*R(2,3), R(1,3)*R(2,1), R(1,1)*R(2,2)];
        
        Qb=[0, 0, 0, R(1,3)*R(1,2), R(1,1)*R(1,3), R(1,2)*R(1,1);
            0, 0, 0, R(2,3)*R(2,2), R(2,1)*R(2,3), R(2,2)*R(2,1);
            0, 0, 0, R(3,3)*R(3,2), R(3,1)*R(3,3), R(3,2)*R(3,1);
            0, 0, 0, R(2,3)*R(3,2), R(2,1)*R(3,3), R(2,2)*R(3,1);
            0, 0, 0, R(3,3)*R(1,2), R(3,1)*R(1,3), R(3,2)*R(1,1);
            0, 0, 0, R(1,3)*R(2,2), R(1,1)*R(2,3), R(1,2)*R(2,1)];
        
        Q  = Qa+Qb;  % rotation matrix for   symmetric tensors
        
%       A4 = Q *A4*R4*inv(Q)*Id4;  % Fourth-order tensor in lab coords.
        A4 = Q * A4 * Q';          % From Nadeau and Ferrari, 1998.


end  % ends the hybrid/orthotropic decision block from the top of the file
