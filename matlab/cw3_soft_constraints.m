load CraneParameters;
xRange = [0 0.52];
yRange = [0 0.62];
Ts=1/10;
Tf=2; % duration of prediction horizon in seconds
N=ceil(Tf/Ts); % ceiling to ensure horizon length N is an integer

%% Load the dynamics matrices using a solution from last assignment
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);

%% Define other simulation parameters
T=4; % duration of simulation
xTarget=0.4*[xRange(2) 0 yRange(2) 0 0 0 0 0]'; % target equilibrium state
xZero=xRange(1); yZero=yRange(1);
x0=[xZero 0 yZero 0 0 0 0 0]'; % initial state

%% Declare penalty matrices and tune them here:
Q=zeros(8);
Q(1,1)=2; % weight on X
Q(3,3)=3; % weight on Y
% In order to test constraints, do not penalise angle or derivative of
% angle. Instead impose soft constraints on them.
%Q(5,5)=1; % weight on theta
%Q(7,7)=1; % weight on psi
R=eye(2)*0.02; % very small penalty on input to demonstrate hard constraints
P=Q; % terminal weight

%% Declare contraints
% Declaring constraints only on states (X,Y,theta,psi) and inputs u
angleConstraint=3*pi/180; % in radians
cl=[0;  0; -angleConstraint;  -angleConstraint];
ch=[0.8*xRange(2);  0.8*yRange(2);  angleConstraint;  angleConstraint  ];
ul=[-1; -1];
uh=[1; 1];
% constrained vector is Dx, hence
D=zeros(4,8);D(1,1)=1;D(2,3)=1;D(3,5)=1;D(4,7)=1;

%% Compute stage constraint matrices and vector
[Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh);

%% Compute trajectory constraints matrices and vector
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);

%% Compute QP constraint matrices
[Gamma,Phi] = genPrediction(A,B,N); % get prediction matrices:
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);

%% Compute QP cost matrices
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);

%% Compute matrices and vectors for soft constraints
% Define weights for constraint violations
rho = 2e3; % weight for exact penalty term
S = 1.4e-3*eye(size(D,1)); % small positive definite quadratic cost to ensure uniqueness
[Hs,gs,Fs,bs,Js,Ls] = genSoftPadding(H,F,bb,J,L,S,rho,size(B,2));