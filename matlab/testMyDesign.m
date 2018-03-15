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
%         c = shape([0.05, 0.05],[0.415,0.415],0.02,0,0);
    case 2
        c = [0.00, 0.05;
             0.25, 0.30;
             0.50, 0.05;
             0.45, 0.00;
             0.25, 0.20;
             0.05, 0.00];
%         c = shape([0.05 0.05], [0.45, 0.05], 0.01, [0.25 0.25], 0);       
%         c = shape([0.05 0.05], [0.45, 0.45], 0.1, [0.05 0.45], 0);       
end


%% Some other parameters
% The starting point
startingPoint = [0.05, 0.05];

% The target point
switch (testShape)
    case 1
        targetPoint = [0.4, 0.4];
    case 2
        targetPoint = [0.45, 0.05];
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


%% Save the data for the simulation 
save('workspace.mat');


%% Open the model
simModel = 'SimscapeCrane_ClosedLoop';
open(simModel); 


%% Import the data into Simulink
mws  = get_param(simModel, 'modelworkspace');
mws.DataSource = 'MAT-File';
mws.FileName = 'workspace';
mws.reload();


%% Setup the simulation time
set_param(bdroot, 'StopTime', num2str(60) );


%% Update the controller blocks
contH = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', [simModel, '/MPController']);
estiH = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', [simModel, '/State_Estimator']);
tgenH = find(slroot, '-isa', 'Stateflow.EMChart', 'Path', [simModel, '/Target_Generator']);

% Error test the block handles
if ( isempty(contH) || isempty(estiH) || isempty(tgenH) )
    error('Unable to get block handle');
end

% Read the controller, state estimator, and target generator into Simulink
contH.Script = fileread('myMPController.m');
estiH.Script = fileread('myStateEstimator.m');
tgenH.Script = fileread('myTargetGenerator.m');


%% Run the actual simulation
sim(simModel);


%% Visualize the Course
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