%% Modify the following function for your setup function
% DONE #############################################################
% IDEA: take c in pieces and work on optimising a route piece by piece.
% Then, for all but the last set of constraints, apply a weak penalty on 
% the final position. This way, the controller will not spend too long
% trying to perfectly fit intermediary positions before moving on.
% ##################################################################
% NEXT STEPS:
% Soft constraints on shape 2. make the constraints the middle path through
% the shape, and then let the soft constraints permit the system to deviate
% slightly.
% TRY SETTING X_HAT TO THE COORDINATES OF THE PENDULUM RATHER THAN THAT OF
% THE CART!!
function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)
    trackwidth = c(2,2) - c(5,2);
    utol = 0.15*trackwidth; ltol = 0.15*trackwidth;
    angleConstraint = 2*pi/180; % in radians
    midpoint = 0.5; % distance between 2 mid points to set 1st target
    % Input constraints (hard)
    inputAttenuation = 0.8;
    ul=inputAttenuation*[-1; -1];
    uh=inputAttenuation*[1; 1];    
    param.x_star = c(2,1) - 0.01;  
    % This is a sample way to send reference points
    % Set targets
    offsetTP1 = [0.0 0.0];
    param.TP1 = midpoint*(c(2,:)+offsetTP1) + (1-midpoint)*(c(5,:)+offsetTP1);
    param.TP2 = targetPoint;
    
    param.rTol = eps_r;
    param.tTol = eps_t;
    
    
    load CraneParameters;
    Ts=1/20;
    Tf=2; % duration of prediction horizon in seconds
    N=ceil(Tf/Ts);
    [A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);    
    %% Declare penalty matrices and tune them here:
%     Q = diag([10 0 10 0 0 0 0 0]);
    Q = diag([1 0 1 0 1 0 1 0]);

    R = eye(2)*0.01; % very small penalty on input to demonstrate hard constraints
    P = Q; % terminal weight

    %% Split shape into 2 quadrilaterals
    c1 = c([1 2 5 6], :);
    c2 = c([2 3 4 5], :);
    % Skew the midpoint so that it doesn't have infinite gradient
    c1(3,1) = c1(3,1) + 0.02;
    c2(4,1) = c2(4,1) - 0.02;
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
    if x_hat(1) < param.x_star
        r(1,1) = param.TP1(1);
        r(3,1) = param.TP1(2);
    else
        condition = (abs(x_hat(1) - param.TP2(1)) < param.tTol) &...
                    (abs(x_hat(3) - param.TP2(2)) < param.tTol) &...
                    (abs(x_hat(2)) < param.rTol) &...
                    (abs(x_hat(4)) < param.rTol)&...
                    (abs(x_hat(6)) < param.rTol)&...
                    (abs(x_hat(8)) < param.rTol);
        if condition
            r(1:8) = x_hat(1:8);
        else
            r(1,1) = param.TP2(1);
            r(3,1) = param.TP2(2);
        end
    end   
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
    x_hat(1:8) = (param.C^-1)*y;
end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% MPC Controller
    opt = mpcqpsolverOptions;
    opt.IntegrityChecks = false;%% for code generation
    opt.FeasibilityTol = 1e-3;
    opt.DataType = 'double';    
    %% Check if we crossed the turning point yet
    if x_hat(1) < param.x_star
        w = x_hat(1:8) - r(1:8);
        f = w'*param.G1';
        b = -(param.bb1 + param.J1*x_hat(1:8) + param.L1*r(1:8));
        [v, ~, ~, ~] = mpcqpsolver(param.H1, f', -param.F1, b, [], zeros(0,1), false(size(param.bb1)), opt);        
    else
        w = x_hat(1:8) - r(1:8);
        f = w'*param.G2';
        b = -(param.bb2 + param.J2*x_hat(1:8) + param.L2*r(1:8));
        [v, ~, ~, ~] = mpcqpsolver(param.H2, f', -param.F2, b, [], zeros(0,1), false(size(param.bb2)), opt);        
    end   
    u = v(1:2);
end % End of myMPController
