clear all
rng default
load CraneParameters;
xRange = [0 0.52];
yRange = [0 0.62];
Ts=1/10;
Tf=2; % duration of prediction horizon in seconds
N=ceil(Tf/Ts);
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

%% Declare penalty matrices:
Q=diag(1:8);
P=diag(3:10);
R=diag([2;3]);

%% Declare contraints
angleConstraint=3*pi/180; % in radians
cl=[0;  0; -angleConstraint;  -angleConstraint];
ch=[0.8*xRange(2);  0.8*yRange(2);  angleConstraint;  angleConstraint];
ul=[-1; -1];
uh=[1; 1];
% constrained vector is Dx, hence
D=zeros(4,8);D(1,1)=1;D(2,3)=1;D(3,5)=1;D(4,7)=1;

%% Compute stage constraint matrices and vector
[Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh);

%% Compute trajectory constraints matrices and vector
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

%% Compute QP constraint matrices
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

%% Compute QP cost matrices
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
% Prepare cost and constraint matrices for mpcqpsolver
% Calculating the inverse of the lower triangular H. see doc mpcqpsolver.
H = chol(H,'lower');
H=(H'\eye(size(H)))';

%% Run a linear simulation to test your genMPController function
xTarget=[0.3 0 0.45 0 0 0 0 0]';% target equilibrium state
x0=[0.52/2 0 0.62/2 0 0 0 0 0]'; % starting offset
T = 5; % simulation duration
iA = false(size(bb));
t=0:Ts:T;
x=[x0, zeros(8,length(t)-1)];
for t_step=1:length(t)-1
    [u,status,iA] = genMPController(H,G,F,bb,J,L,x(:,t_step),xTarget,size(B,2),iA);
    x(:,t_step+1)=A*x(:,t_step)+B*u;
    if status ~= -1
    %   status = -1 --> Cannot solve QP subject to constraints provided     
        disp(status);
    end
end

%% Plot results
allTitles=[{'X'},{'dX/dt'},{'Y'},{'dY/dt'},{'\theta'},{'d\theta/dt'},{'\psi'},{'d\psi/dt'}];
figure;
for i=1:8
     subplot(4,2,i);
     plot(t,x(i,:));
     title(allTitles(i));
end




