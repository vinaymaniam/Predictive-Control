function r = myTargetGenerator(x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % Make the crane go to (xTar, yTar)
    r(1,1) = param.TP(1);
    r(3,1) = param.TP(2);
%     if param.useDistRej
%         r = param.DR1\(param.DR2*[x_hat(9:16); r(1:8)]);
%     end
    
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
