%% Modify the following function for your setup function
% IDEA: take c in pieces and work on optimising a route piece by piece.
% Then, for all but the last set of constraints, apply a weak penalty on 
% the final position. This way, the controller will not spend too long
% trying to perfectly fit intermediary positions before moving on.


% NEED TO PASS PARAM BY REFERENCE IN MYMPCCONTROLLER ELSE IT WONT UPDATE
function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
    % This is a sample static K matrix for the controller
    param.K = [1, 0, 0, 0, 0, 0, 0, 0;
               0, 0, 1, 0, 0, 0, 0, 0];
       
    % This is a sample way to send reference points
    param.OTP = targetPoint(1:2);
    
    param.er = eps_r;
    param.et = eps_t;
    
    param.start = startingPoint;
    if size(c,1) == 6
        param.c = c(2:5,:);
        param.c(4,1) = param.c(4,1) - 0.02;
        c = c([1 2 5 6],:);
        c(3,1) = c(3,1) + 0.03;
        param.TP = 0.5*(c(2,:) + c(3,:));
        param.mode = 0;
    else
        param.TP = param.OTP;
        param.mode = 1;        
    end    
    
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
    
    tol = 0.00;
    
    for i = 1:size(c,1)
        i2 = mod(i+1,size(c,1));
        if i2 == 0
            i2 = size(c,1);
        end        
        if c(i,1) > c(i2,1)
            modifier = -1;
            c1y = c(i,2) + tol;
            c2y = c(i2,2) + tol;
        else
            modifier = 1;
            c1y = c(i,2) - tol;
            c2y = c(i2,2) - tol;            
        end
%         coeff = polyfit([c(i,1), c(i2,1)], [c(i,2), c(i2,2)], 1);
        coeff = polyfit([c(i,1), c(i2,1)], [c1y, c2y], 1);
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
%     param.gs = gs;
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
    
        
    param.x_star = c(2,1) - 0.02;       
end % End of mySetup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your target generation
function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % This is the block in which I have full control over how i decide where the crane
    % should go.
    % Make the crane go to (xTar, yTar)
    r(1,1) = param.TP(1);
    r(3,1) = param.TP(2);
    r = r(1:8);
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %%    
%     % By default, just pass the system measurements through
%     x_hat( 1:length(y),1 ) = y;
    x_hat = (param.C^-1)*y;
end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% Check if we crossed the turning point yet
    if (param.mode == 0) & (x_hat(1) > param.x_star)
        param = mySetup(param.c, x_hat([]), param.OTP, param.er, param.et);
    end
    %% MPC Controller
    opt = mpcqpsolverOptions;
    opt.IntegrityChecks = false;%% for code generation
    opt.FeasibilityTol = 1e-3;
    opt.DataType = 'double';
    %% your code starts here
    % Cholksey and inverse already computed and stored in H
    w = x_hat - r;
    f = w'*param.G';

    b = -(param.bb + param.J*x_hat + param.L*r);
    [v, ~, ~, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), false(size(param.bb)), opt);
%     w = x_hat - r;
%     f = [w'*param.G', param.gs'];
% 
%     b = -(param.bb + param.J*x_hat + param.L*r);
%     [v, ~, ~, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), false(size(param.bb)), opt);
    %% your remaining code here
    u = v(1:2);
end % End of myMPController