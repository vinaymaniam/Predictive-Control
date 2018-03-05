function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %% Pendulum is assumed to be of length 0.47m
    x_hat(1:8) = param.C\y;
end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
