function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % This is the block in which I have full control over how i decide where the crane
    % should go.
    % Make the crane go to (xTar, yTar)
    if x_hat(1) < param.x_star
        r(1,1) = param.TP1(1);
        r(3,1) = param.TP1(2);
    else
        r(1,1) = param.TP2(1);
        r(3,1) = param.TP2(2);
    end
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
