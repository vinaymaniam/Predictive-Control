function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)

    % This is a sample static K matrix for the controller
    param.K = [1, 0, 0, 0, 0, 0, 0, 0;
               0, 0, 1, 0, 0, 0, 0, 0];

    % This is a sample way to send reference points
    param.xTar = targetPoint(1);
    param.yTar = targetPoint(2);
    
    param.rTol = eps_r;
    param.tTol = eps_t;
    
    param.start = startingPoint;
    
    load CraneParameters;
    Ts=1/10;
    Tf=2; % duration of prediction horizon in seconds
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);    
    %% Declare penalty matrices and tune them here:
    Q=zeros(8);
    Q(1,1)=10; % weight on X
    Q(3,3)=10; % weight on Y

    R=eye(2)*0.01; % very small penalty on input to demonstrate hard constraints
    P=Q; % terminal weight

    %% Declare contraints
    % Constrain only states (X,Y,theta,psi)
    % Constrained vector is Dx, hence
    %% Construct constraint matrix D
    % General form
    D = zeros(size(c,1) + 2, 8);
    ch = zeros(size(c,1) + 2, 1);
    for i = 1:size(c,1)
        i2 = mod(i+1,size(c,1));
        if i2 == 0
            i2 = size(c,1);
        end
        coeff = polyfit([c(i,1), c(i2,1)], [c(i,2), c(i2,2)], 1);
        if c(i,1) > c(i2,1)
            modifier = -1;
        else
            modifier = 1;
        end
        D(i,1) = coeff(1) * modifier*(-1);
        D(i,3) = modifier;     
        ch(i) = modifier * coeff(2);
    end
    angleConstraint=2*pi/180; % in radians
    D(end-1,5) = 1;      ch(end-1) = angleConstraint;
    D(end,7) = 1;        ch(end) = angleConstraint;  

    %% End of construction        
    % Input constraints (hard)
    ul=[-1; -1];
    uh=[1; 1];

    %% Compute stage constraint matrices and vector
%     [Dt,Et,bt]=genStageConstraints(A,B,D,[],constraints,ul,uh);
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

%     %% Compute matrices and vectors for soft constraints
%     % Define weights for constraint violations
%     rho = 1e3; % weight for exact penalty term
%     S = 1e-3*eye(size(D,1)); % small positive definite quadratic cost to ensure uniqueness
%     [Hs,gs,Fs,bs,Js,Ls] = genSoftPadding(H,F,bb,J,L,S,rho,size(B,2));
% 
%     %% replace matrices and vectors to simplify code
%     H = Hs; F = Fs; bb = bs; J = Js; L = Ls;

    %% Prepare cost and constraint matrices for mpcqpsolver
    H = chol(H,'lower');
    H=H\eye(size(H));    
    
    param.H = H;
    param.G = G;
%     param.gs = gs;
    param.F = F;
    param.J = J;
    param.L = L;
    param.bb = bb;
    
    param.A = A;
    param.C = C;
    
end % End of mySetup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your target generation