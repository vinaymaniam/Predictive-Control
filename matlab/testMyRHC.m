clear variables
close all
%% Load the parameters
load('Params_Simscape.mat');
load('SSmodelParams.mat');
%% Declare simulation parameters
Ts=1/20;
N=ceil(3/Ts);
T=20;
xTarget=[0.4 0 0.5 0 0 0 0 0]';
% Load the matrices using a solution from the previous assignment
[A,B,C,D] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);
%% Declare penalty matrices and tune them here:
Q=eye(8);
R=eye(2);
P=Q;  
%% Compose prediction matrices for RHC
% your myPrediction function is called here and the variables Gamma and Phi are 
% declared in the workspace. 
[Gamma,Phi]=genPrediction(A,B,N);
%% Declare RHC control law
% The linear control law K is declared here. It will be visible to a Simulink 
% constant block and can be used to implement the control law.
% See how this is implemented here: SimscapeCrane_RHC/Controllers
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);
K = genRHC(H,G,size(B,2));
%% Run the simulations for your controller and the PID controller
% Select controller, Uncomment as required
% controlCase=1; % your RHC 
% controlCase=2; % P(ID)controller

% Open the model
SimscapeCrane_RHC; 

controlCase=2;
sim('SimscapeCrane_RHC');
responsePP.output=GantryCraneOutput;
responsePP.input=GantryCraneInput;

controlCase=1;
sim('SimscapeCrane_RHC');
responseRHC.output=GantryCraneOutput;
responseRHC.input=GantryCraneInput;

%% visualise the performance:
help GantryResponsePlot
GantryResponsePlot(responsePP.output.time,responsePP.input.signals.values,...
    responsePP.output.signals.values,[-1 -1],[1 1],[0 0],[xRange(2) yRange(2)],[1 3],xTarget,'PID performance');
GantryResponsePlot(responseRHC.output.time,responseRHC.input.signals.values,...
    responseRHC.output.signals.values,[-1 -1],[1 1],[0 0],[xRange(2) yRange(2)],[1 3],xTarget,'RHC performance');