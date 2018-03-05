function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
    trackwidth = c(2,2) - c(5,2);
    utol = 0.2*trackwidth; ltol = 0.2*trackwidth;
    angleConstraint = 2*pi/180; % in radians
    midpoint = 0.3; % distance between 2 mid points to set 1st target
    % Input constraints (hard)
    inputAttenuation = 1;
    ul=inputAttenuation*[-1; -1];
    uh=inputAttenuation*[1; 1];     
    % This is a sample way to send reference points
    % Set targets
%     param.TP1 = 0.5*c(2,:) + 0.5*c(5,:);
    param.TP1 = midpoint*c(2,:) + (1-midpoint)*c(5,:);
    param.TP2 = targetPoint;
    
    param.er = eps_r;
    param.et = eps_t;
    
    
    load CraneParameters;
    Ts=1/20;
    Tf=2; % duration of prediction horizon in seconds
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);    
    %% Declare penalty matrices and tune them here:
%     Q = diag([10 0 10 0 1 0 1 0]);
    Q = diag([30 0.001 30 0.001 1 0.001 1 0.001]);

    R = eye(2)*0.5; % very small penalty on input to demonstrate hard constraints
    P = Q; % terminal weight

    %% Split shape into 2 quadrilaterals
    c1 = c([1 2 5 6], :);
    c2 = c([2 3 4 5], :);
    % Skew the midpoint so that it doesn't have infinite gradient
    c1(3,1) = c1(3,1) + 0.01;
    c2(4,1) = c2(4,1) - 0.01;
%%  c1
    %% Construct constraint matrix D
    % General form
    D = zeros(size(c1,1) + 2, 8);
    ch = zeros(size(c1,1) + 2, 1);        
    
    for i = 1:size(c1,1)
        i2 = mod(i+1,size(c1,1));
        if i2 == 0
            i2 = size(c1,1);
        end        
        if c1(i,1) > c1(i2,1)
            modifier = -1;
            c1y = c1(i,2) + ltol;
            c2y = c1(i2,2) + ltol;
        else
            modifier = 1;
            c1y = c1(i,2) - utol;
            c2y = c1(i2,2) - utol;            
        end
        coeff = polyfit([c1(i,1), c1(i2,1)], [c1y, c2y], 1);
        D(i,1) = coeff(1) * modifier*(-1);
        D(i,3) = modifier;     
        ch(i) = modifier * coeff(2);
    end    
    D(end-1,5) = 1;      ch(end-1) = angleConstraint;
    D(end,7) = 1;        ch(end) = angleConstraint;  

    %% End of construction        

    %% Compute stage constraint matrices and vector
    DA = D*A;
    DB = D*B;
    I = eye(size(B,2));
    O = zeros(size(B,2),size(A,2));
    Dt = [DA; O; O];
    Et = [DB; I; -I];
    bt = [ch; uh; -ul];    
    %% Compute trajectory constraints matrices and vector
    [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

    %% Compute QP constraint matrices
    [Gamma,Phi] = genPrediction(A,B,N); % get prediction matrices:
    [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

    %% Compute QP cost matrices
    [H,G]=genCostMatrices(Gamma,Phi,Q,R,P,N);
    
    %% Prepare cost and constraint matrices for mpcqpsolver
    H = chol(H,'lower');
    H=H\eye(size(H));    
    
    param.H1 = H;
    param.G1 = G;    
    param.F1 = F;
    param.J1 = J;
    param.L1 = L;
    param.bb1 = bb;
%%  c2
    %% Construct constraint matrix D
    % General form
    D = zeros(size(c2,1) + 2, 8);
    ch = zeros(size(c2,1) + 2, 1);        
    
    for i = 1:size(c2,1)
        i2 = mod(i+1,size(c2,1));
        if i2 == 0
            i2 = size(c2,1);
        end        
        if c2(i,1) > c2(i2,1)
            modifier = -1;
            c1y = c2(i,2) + ltol;
            c2y = c2(i2,2) + ltol;
        else
            modifier = 1;
            c1y = c2(i,2) - utol;
            c2y = c2(i2,2) - utol;            
        end
        coeff = polyfit([c2(i,1), c2(i2,1)], [c1y, c2y], 1);
        D(i,1) = coeff(1) * modifier*(-1);
        D(i,3) = modifier;     
        ch(i) = modifier * coeff(2);
    end
    angleConstraint = 2*pi/180; % in radians
    D(end-1,5) = 1;      ch(end-1) = angleConstraint;
    D(end,7) = 1;        ch(end) = angleConstraint;  

    %% End of construction        
    % Input constraints (hard)
    ul=[-1; -1];
    uh=[1; 1];

    %% Compute stage constraint matrices and vector
    DA = D*A;
    DB = D*B;
    I = eye(size(B,2));
    O = zeros(size(B,2),size(A,2));
    Dt = [DA; O; O];
    Et = [DB; I; -I];
    bt = [ch; uh; -ul];    
    %% Compute trajectory constraints matrices and vector
    [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

    %% Compute QP constraint matrices
    [Gamma,Phi] = genPrediction(A,B,N); % get prediction matrices:
    [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

    %% Compute QP cost matrices
    [H,G]=genCostMatrices(Gamma,Phi,Q,R,P,N);
    
    %% Prepare cost and constraint matrices for mpcqpsolver
    H = chol(H,'lower');
    H=H\eye(size(H));    
    
    param.H2 = H;
    param.G2 = G;    
    param.F2 = F;
    param.J2 = J;
    param.L2 = L;
    param.bb2 = bb;
%%  Rest of code
    param.A = A;
    param.C = C;
    
        
    param.x_star = c(2,1) - 0;       
end % End of mySetup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your target generation
