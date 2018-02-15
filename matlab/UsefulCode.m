% This script serves as a collection of sample code to achieve some
% tasks you  might encounter throughout the course. 
%% Plot step input response of a system
% generate a random system
randomSystem=rss(5,3,2);
% sample the system at 1/25
discRandSystem=c2d(randomSystem,1/25);
% plot step response
step(discRandSystem);

%% Compute output of a discrete-time model provided an input vector
in=rand(100,1);
Ts=1/10;
randomSystem=rss(5,3,1);
discRandSystem=c2d(randomSystem,Ts);
lsimplot(discRandSystem,in,Ts:Ts:10);
%% Play around with the simulink model:
% declare simulation parameters
Ts=1/10;
T=20;
% Load physical modelling parameters
load('Params_Simscape.mat');
%% Using data inspector:
% http://www.mathworks.com/help/simulink/ug/simulation-data-inspector-overview.html