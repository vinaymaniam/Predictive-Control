%% Modify the following function for your setup function
%% 3.15 seconds!!!!!!
function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%%  Configure here    
    trackwidth = sqrt(sum((c(2,:) - c(3,:)).^2));    
    tol = 0.15*trackwidth;    
    inputAttenuation = 1;%0.78; % best=0.78
    ul=inputAttenuation*[-1; -1];
    uh=inputAttenuation*[1; 1];    
    %% Choose which modifications to use
    param.mod = 0; % 0 = NO MOD, 1 = OFFSET BLOCKING
    param.soft = 0; % 0 for hard, 1 for soft
    useRatePen = 1; % 0 for no rate penalties to input    
    param.useDistRej = 0; % 0 for no disturbance rejection   
    %%        
    angleConstraint = 4*pi/180; % in radians
%   END OF CONFIGURATION    
    param.TP = targetPoint;
    
    param.rTol = eps_r;
    param.tTol = eps_t;
    
    param.start = startingPoint;
    
%     load CraneParameters;
    load SSmodelParams.mat;
    load Params_Simscape.mat;
    Ts=1/20;    
    Tf=2; % duration of prediction horizon in seconds
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);  
    param.A = A;
    param.B = B;
    param.C = C;
    
    %% Declare penalty matrices and tune them here:
    Q=zeros(8);
    % penalties = [10,0,10,0,50,0,50,0]; % works pretty well
    pos = 3; vel = 0; angl = 20; rangl = 0;%0.03 no I/P rate pen 
%     pos = 10; vel = 0; angl = 50; rangl = 0;%0.03 no I/P rate pen 
    penalties = [pos,vel,pos,vel,angl,rangl,angl,rangl];
    for i = 1:length(penalties)
        Q(i,i) = penalties(i);
    end
    %% CHANGE TO VERY LOW NUMBER(EG. 0.0001) AND ADD RATE PENALTIES
    R=eye(2)*0.0001; % very small penalty on input to demonstrate hard constraints
%     R=eye(2)*0.003; % very small penalty on input to demonstrate hard constraints
    P=Q; % terminal weight
    %% Smart Choice of P
    [K,~,~] = dlqr(A, B, Q, R);
    P = dlyap((A-B*K)', Q + K'*R*K);    
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
    D = zeros(size(c,1) + 2, 8);
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
    
    %% Rate penalties
    if useRatePen == 1
        R2 = 0.0002*eye(2);
        T = [zeros(2,(N-1)*2), zeros(2,2);
             eye((N-1)*2),     zeros((N-1)*2,2)];
        RatePenMat = ((eye(N*2)) - T)'*kron(eye(N),R2)*((eye(N*2)) - T);
        H = H + 2*RatePenMat;
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
    
    %% Disturbance Rejection stuff
    Cd = 0.01*eye(8);
    Bd = 0.01*eye(8);
    Hdr = eye(8);
    Mdrd = Hdr * Cd;
    Mdr = Hdr * C;
    param.DR1 = [eye(size(A,1))-A,                   -param.B;
                 Mdr,      zeros(size(Mdr,1),size(param.B,2))];
    param.DR2 = [Bd,      zeros(size(Bd));
                 -Mdrd,   eye(size(Mdrd,1))];   
    param.Adr = [A,                Bd;
                 zeros(size(A,1)), eye(size(A,1))];
    param.Bdr = [B; zeros(8,2)];
    
    Ct = [param.C, Cd]';
    Qt = 0.1*eye(16);
    Rv = 0.1*eye(8); 
    [Sigma,~,~] = dare(param.Adr',Ct,Qt,Rv);
    param.Ldr = param.Adr'*Sigma*Ct/((Ct'*Sigma*Ct + Rv));
    param.Cd = Cd;
    param.startingPoint = startingPoint;
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
    if param.useDistRej
        r = param.DR1\(param.DR2*[x_hat(9:16); r(1:8)]);
    end
    
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %% Pendulum is assumed to be of length 0.47m
    x_hat(1:8) = param.C\y;    
    %% Disturbance Rejection
    persistent state;     
    if isempty(state)
        state = zeros(16,1);
        state(1) = param.startingPoint(1);
        state(3) = param.startingPoint(2);
    else
        if param.useDistRej == 1                                          
            state = param.Adr * state + param.Bdr*u +...
                param.Ldr*(y - [param.C, param.Cd]*state);
        end
    end
    if param.useDistRej == 1
        x_hat = state;
    end
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
    persistent iA;
    if isempty(iA)
        iA = false(size(param.bb));
    end
    if ~condition    
        %% MPC Controller
        opt = mpcqpsolverOptions;
        opt.IntegrityChecks = false;%% for code generation
        opt.FeasibilityTol = 1e-3;
        opt.DataType = 'double';
        %% your code starts here            
        % Cholksey and inverse already computed and stored in H
        w = x_hat(1:8) - r(1:8);
        
        f = w'*param.G'; % Hard
        b = -(param.bb + param.J*x_hat(1:8) + param.L*r(1:8));
        [ubar, ~, iA, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), iA, opt);
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
