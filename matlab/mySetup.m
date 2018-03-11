function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
%%  Configure here
    param.mod = 0; % 0 = NO MOD, 1 = OFFSET BLOCKING
    tol = 0.01;    
    inputAttenuation = 0.78; % best=0.78
    ul=inputAttenuation*[-1; -1];
    uh=inputAttenuation*[1; 1];
    Tf=1.5; % duration of prediction horizon in seconds
    param.soft = 0; % 0 for hard, 1 for soft
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
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);    
    %% Declare penalty matrices and tune them here:
    Q=zeros(8);
    Q = diag([10,0,10,0,50,0,50,0]);
    R=eye(2)*0.01; % very small penalty on input to demonstrate hard constraints
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
%     % Remove redundant lines from D <--- didnt really make it faster
%     idx = lblines(round(lblines(:,1),3) == round(ublines(:,1),3),3);
%     D2 = D; ch2 = ch; cl2 = cl;
%     D = []; ch = []; cl = [];
%     for i=1:size(D2,1)
%         disp(any(idx == i))
%         if ~any(idx == i)
%             D = [D; D2(i,:)];
%             ch = [ch; ch2(i,:)];
%             cl = [cl; cl2(i,:)];
%         else
%             cl = [cl; -lblines(lblines(:,3)==i,2)];
%         end            
%     end
    
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
