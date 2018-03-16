function r = myTargetGenerator(x_hat, param)

    %% Do not delete this line
    % Create the output array of the appropriate size
    r = zeros(10,1);
    %%
    % This is the block in which I have full control over how i decide where the crane
    % should go.
    % Make the crane go to (xTar, yTar)
    persistent stuck;
    persistent stuck2;
    if isempty(stuck)
        stuck = 0;
        stuck2 = 0;
    end
    if (param.toggle*x_hat(3) < param.toggle*(x_hat(1)*param.switch_line(1) + param.switch_line(2))) && stuck2 == 0
        r(1,1) = param.TP1(1);
        r(3,1) = param.TP1(2);
    else        
        r(1,1) = param.TP2(1);
        r(3,1) = param.TP2(2);
    end    
    radius = sqrt((abs(x_hat(1) - param.TP1(1)))^2+(abs(x_hat(3) - param.TP1(2)))^2);
    % 0.0003 was manually selected from data collection of the oscillatory
    % response of the overly aggressive system around TP
    if (radius < 0.003) && (stuck2 == 0)
        stuck = stuck + 1;        
        if stuck > 2
            r(1,1) = param.TP2(1);
            r(3,1) = param.TP2(2);            
            stuck2 = 1;
        end        
    end
end % End of myTargetGenerator


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Modify the following function for your state estimator (if desired)
