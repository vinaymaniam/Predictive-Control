%% Modify the following function for your setup function
function [ param ] = mySetup(c, startingPoint, targetPoint, eps_r, eps_t)

    % This is a sample static K matrix for the controller
    param.K = [1, 0, 0, 0, 0, 0, 0, 0;
               0, 0, 1, 0, 0, 0, 0, 0];

    % This is a sample way to send reference points
    param.xTar = targetPoint(1);
    param.yTar = targetPoint(2);

end % End of mySetup


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your target generation
function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%

    % Make the crane go to (xTar, yTar)
    r(1,1) = param.xTar;
    r(3,1) = param.yTar;

end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %%

    % By default, just pass the system measurements through
    x_hat( 1:length(y),1 ) = y;

end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
function u = myMPController(r, x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %%

    % Do a stateic gain state feedback controller
    u = param.K*(r(1:8,1) - x_hat(1:8,1));

end % End of myMPController
