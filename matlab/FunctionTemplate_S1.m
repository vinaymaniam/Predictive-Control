%% Modify the following function for your setup function

function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%%  Configure here
    param.mod = 0; % 0 = NO MOD, 1 = OFFSET BLOCKING
    tol = 0.01;    
    inputAttenuation = 1;%0.78; % best=0.78
    ul=inputAttenuation*[-1; -1];
    uh=inputAttenuation*[1; 1];    
    param.soft = 0; % 0 for hard, 1 for soft
    
    useDistRej = 0; % 0 for no disturbance rejection    
    
    if param.soft == 0
        angleConstraint = 8*pi/180; % in radians
    else
        angleConstraint = 2*pi/180; % in radians
    end
%   END OF CONFIGURATION    
    param.TP = targetPoint;
    
    param.rTol = eps_r;
    param.tTol = eps_t;
    
    param.start = startingPoint;
    
    load CraneParameters;
    Ts=1/20;    
    Tf=1.5; % duration of prediction horizon in seconds
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);  
    %% Disturbance Rejection
    if useDistRej
        nd = 8; % number of features used to model disturbance
        Bd = [eye(nd); zeros(size(A,1)-nd,nd)];
        A = [A,                             Bd; 
             zeros(size(Bd,2), size(A,2)),  eye(size(Bd,2))];
        B = [B; zeros(size(Bd,2),size(B,2))];
        Cd = [eye(nd); zeros(size(C,1)-nd, nd)];
        C = [C, Cd];
        xsz = 8 + nd;
    else
        xsz = 8;
    end
    
    %% Declare penalty matrices and tune them here:
    Q=zeros(xsz);
    penalties = [10,0,10,0,50,0,50,0];
    for i = 1:length(penalties)
        Q(i,i) = penalties(i);
    end
    R=eye(2)*0.003; % very small penalty on input to demonstrate hard constraints
    P=Q; % terminal weight
    if param.mod == 1
        % Ricatti Thing    
        [P,~,~] = dare(A,B,Q,R);
        K_lqr = -((B'*P*B + R)^-1)*B'*P*A;
        param.K_lqr = K_lqr;
    else
        param.K_lqr = zeros(2,8);
    end
    %% Construct constraint matrix D
    % General form
    D = zeros(size(c,1) + 2, xsz);
    ch = zeros(size(c,1) + 2, 1);       
    lblines = [];
    ublines = [];
    for i = 1:size(c,1)
        i2 = mod(i+1,size(c,1));
        if i2 == 0
            i2 = size(c,1);
        end        
        if c(i,1) > c(i2,1)
            modifier = -1;
            c1y = c(i,2) + tol;
            c2y = c(i2,2) + tol;
            coeff = polyfit([c(i,1), c(i2,1)], [c1y, c2y], 1);
            lblines = [lblines; [coeff, i]];
        else
            modifier = 1;
            c1y = c(i,2) - tol;
            c2y = c(i2,2) - tol;    
            coeff = polyfit([c(i,1), c(i2,1)], [c1y, c2y], 1);
            ublines = [ublines; [coeff, i]];
        end
        D(i,1) = coeff(1) * modifier*(-1);
        D(i,3) = modifier;     
        ch(i) = modifier * coeff(2);
    end    
    D(end-1,5) = 1;      D(end,7) = 1;        
    ch(end-1) = angleConstraint;    
    ch(end) = angleConstraint;      
    cl = -inf*ones(size(ch));
    cl(end-1) = -angleConstraint;
    cl(end) = -angleConstraint;
    
    %% End of construction        
    

    %% Compute stage constraint matrices and vector    
    [Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh);
    %% Compute trajectory constraints matrices and vector
    [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

    %% Compute QP constraint matrices
    [Gamma,Phi] = genPrediction(A,B,N); % get prediction matrices:
    [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

    %% Compute QP cost matrices
    [H,G]=genCostMatrices(Gamma,Phi,Q,R,P,N);
    
    if param.soft == 1
        %% Compute matrices and vectors for soft constraints
        % Define weights for constraint violations
        rho = 1e3; % weight for exact penalty term
        S = 1e-3*eye(size(D,1)); % small positive definite quadratic cost to ensure uniqueness
        [Hs,gs,Fs,bs,Js,Ls] = genSoftPadding(H,F,bb,J,L,S,rho,size(B,2));
        %% replace matrices and vectors to simplify code
        H = Hs; F = Fs; bb = bs; J = Js; L = Ls;
        param.gs = gs;
        param.S = S;
    else
        param.gs = 0;
        param.S = 0;
    end
    %% Offset blocking
    if param.mod == 1
        Iu = eye(N*2); % N*m
        IkronK = kron(eye(N), K_lqr);
        M1 = (Iu - IkronK*Gamma)^-1;
        M2 = IkronK*Phi;
        H = M1'*H*M1;
        G = (G'*M1 + M2'*M1'*H*M1)';
        param.M1 = M1;
        param.M2 = M2;
    else
        param.M1 = 0;
        param.M2 = 0;
    end
    
    %% Prepare cost and constraint matrices for mpcqpsolver
    H = chol(H,'lower');
    H=H\eye(size(H));    
    
    param.H = H;
    param.G = G;    
    param.F = F;
    param.J = J;
    param.L = L;
    param.bb = bb;
    
    param.A = A;
    param.C = C;            
end % End of mySetup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your target generation
function r = myTargetGenerator(x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % Make the crane go to (xTar, yTar)
    r(1,1) = param.TP(1);
    r(3,1) = param.TP(2);
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %% Pendulum is assumed to be of length 0.47m
    x_hat(1:8) = param.C\y;
end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% Check if crane is at target point
    radius = sqrt((abs(x_hat(1) - param.TP(1)))^2+(abs(x_hat(3) - param.TP(2)))^2);
    condition = (radius < param.tTol) &...
                (abs(x_hat(2)) < param.rTol) &...
                (abs(x_hat(4)) < param.rTol)&...
                (abs(x_hat(6)) < param.rTol)&...
                (abs(x_hat(8)) < param.rTol);
%     condition = 0; % Leave the condition to be handled by target gen
    if ~condition    
        %% MPC Controller
        opt = mpcqpsolverOptions;
        opt.IntegrityChecks = false;%% for code generation
        opt.FeasibilityTol = 1e-3;
        opt.DataType = 'double';
        %% your code starts here
        % Cholksey and inverse already computed and stored in H
        w = x_hat(1:8) - r(1:8);
        if param.soft == 0
            f = w'*param.G'; % Hard
        else            
            f = [w'*param.G', param.gs']; % Soft
        end
        b = -(param.bb + param.J*x_hat(1:8) + param.L*r(1:8));
        [ubar, ~, ~, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), false(size(param.bb)), opt);
        %% your remaining code here
        if param.mod == 1
%             ubar = param.M1*param.M2*x_hat(1:8) + param.M1*ubar;
            v = param.K_lqr*x_hat(1:8) + ubar(1:2);
        else
            v = ubar(1:2);
        end
        u = v;
    end
end % End of myMPController
