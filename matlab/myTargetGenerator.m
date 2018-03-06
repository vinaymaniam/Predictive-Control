function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % Make the crane go to (xTar, yTar)
    r(1,1) = param.xTar;
    r(3,1) = param.yTar;
%     condition = (abs(x_hat(1) - param.xTar) < param.tTol) &...
%                 (abs(x_hat(3) - param.yTar) < param.tTol) &...
%                 (abs(x_hat(2)) < param.rTol) &...
%                 (abs(x_hat(4)) < param.rTol)&...
%                 (abs(x_hat(6)) < param.rTol)&...
%                 (abs(x_hat(8)) < param.rTol);
%     if condition
%         r(1:8) = x_hat(1:8);
%     end
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
