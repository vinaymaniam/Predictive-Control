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
