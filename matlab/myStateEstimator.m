function x_hat = myStateEstimator(u, y, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    x_hat = zeros(16,1);
    %% Pendulum is assumed to be of length 0.47m
    x_hat(1:8) = param.C\y;    
    %% Disturbance Rejection
%     persistent state;          
%     if isempty(state)
%         state = zeros(16,1);
%         state(1) = param.startingPoint(1);
%         state(3) = param.startingPoint(2);
%     else
%         if param.useDistRej == 1 
%             x_hat = state;
%             state = param.Adr * state + param.Bdr*u +...
%                 param.Ldr*(y - [param.C, param.Cd]*state);
%         end
%     end
end % End of myStateEstimator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your controller
