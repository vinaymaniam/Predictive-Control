function u = myMPController(r, x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %%

    % Do a stateic gain state feedback controller
    u = param.K*(r(1:8,1) - x_hat(1:8,1));

end % End of myMPController
