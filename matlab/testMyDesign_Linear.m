clear variables
close all

testShape = 1;


%% Create the test shape
switch ( testShape )
    case 1
        c = [0.00, 0.05;
             0.45, 0.50;
             0.50, 0.45;
             0.05, 0.00];
%         c = [0.00, 0.45;
%              0.45, 0.45;
%              0.45, 0.35;
%              0.00, 0.35];
    case 2
        c = [0.00, 0.05;
             0.25, 0.30;
             0.50, 0.05;
             0.45, 0.00;
             0.25, 0.20
             0.05, 0.00];
%         c = [0.00, 0.00;
%              0.00, 0.45;
%              0.45, 0.45;
%              0.45, 0.35;
%              0.10, 0.35;
%              0.10, 0.00];  
end


%% Some other parameters
% The starting point
startingPoint = [0.05, 0.05];
%startingPoint = [0.05, 0.4];

% The target point
switch (testShape)
    case 1
        targetPoint = [0.4, 0.4];
    case 2
        targetPoint = [0.45, 0.05];
        %targetPoint = [0.4, 0.4];
end

% The tolerances
eps_r = 0.02;
eps_t = 0.02;


%% Load the parameters for the Simscape
load('Params_Simscape.mat');
load('SSmodelParams.mat');


%% Extract the student functions
extractFunctions(['FunctionTemplate_S' num2str(testShape), '.m'], 1);


%% Declare other simulation parameters
f = 20;
Ts = 1/f;

xZero = startingPoint(1,1);
yZero = startingPoint(1,2);


%% Call the setup function for the student
param = mySetup(c, startingPoint, targetPoint, eps_r, eps_t);


%% Create the dynamics matrices
[A,B,C,~] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts);
sysd=ss(A,B,C,0,Ts);


%% Set the simulation time
T = 60;


%% Run the actual linear simulation
% create waiting bar
hw=waitbar(0,'Please wait...');
warning('on');

% Initial conditions
x=[xZero; 0; yZero; 0; 0; 0; 0; 0];
y = x;
u = [0; 0];
u_all = [];

% Variable to hold optimization time
numSamples = T/Ts;
allContTime = [];

% Iterate for the simulation time
ctr = 1;
for t=0:Ts:T
    waitbar(t/T,hw,'Please wait...');
    tic
    % Call the state estimator
    x_hat = myStateEstimator(u, y, param);
    
    % Call the target generator
    ref = myTargetGenerator(x_hat, param);
    
    % Call the controller function
    u = myMPController(ref, x_hat, param);    
    contTime=toc;
%     if max(abs(u)) > 0.001
%         disp([u' ctr])
%         ctr = ctr+1;
%     end
    % Simulate
    [y, tt, xx] = lsim(sysd, [u';0 0], [0 Ts], x(:,end));

    % Keep the state variables around for analysis
    x=[x xx(end,:)'];
    u_all = [u_all u];
    
    % Save only the most recent output
    y = y(end,:)';
    
    % Save the computation time for analysis
    allContTime=[allContTime; contTime];
end

close(hw);
t=0:Ts:t;
x=x(:,1:length(t))';


%% Analyze the controller runtime
figure('Name','Controller Runtime');
plot(t, allContTime);
xlabel('Simulation time [s]')
ylabel('Controller runtime [s]')


%% Visualize the Course
GantryCraneOutput.time = t;
GantryCraneOutput.signals.values = x;
GantryCraneInput.signals.values = u_all';

analyzeCourse( GantryCraneOutput, testShape, c, r, startingPoint, targetPoint );
h = circle(targetPoint(1),targetPoint(2),eps_t);
plot(h(:,1),h(:,2),'k');
%% Visualize the performance
ul=[-1; -1];
uh=[1; 1];
cl=[0; 0];
ch=[xRange(2); yRange(2)];
xTarget = [targetPoint(1); 0; targetPoint(2); 0; 0; 0; 0; 0];
GantryResponsePlot(GantryCraneOutput.time,GantryCraneInput.signals.values,...
    GantryCraneOutput.signals.values,ul,uh,cl,ch,[1 3],xTarget,'Nonlinear simulation',eps_r,eps_t);