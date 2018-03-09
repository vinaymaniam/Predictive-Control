function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % This is the block in which I have full control over how i decide where the crane
    % should go.
    % Make the crane go to (xTar, yTar)
    if param.toggle*x_hat(3) < param.toggle*(x_hat(1)*param.switch_line(1) + param.switch_line(2))        
        r(1,1) = param.TP1(1);
        r(3,1) = param.TP1(2);
    else
        condition = (abs(x_hat(1) - param.TP2(1)) < param.tTol) &...
                    (abs(x_hat(3) - param.TP2(2)) < param.tTol) &...
                    (abs(x_hat(2)) < param.rTol) &...
                    (abs(x_hat(4)) < param.rTol)&...
                    (abs(x_hat(6)) < param.rTol)&...
                    (abs(x_hat(8)) < param.rTol);
        if condition
            r(1:8) = x_hat(1:8);
        else
            r(1,1) = param.TP2(1);
            r(3,1) = param.TP2(2);
        end
    end   
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
