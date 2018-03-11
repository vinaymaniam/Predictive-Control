function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% Check if crane is at target point
    radius = sqrt((abs(x_hat(1) - param.TP(1)))^2+(abs(x_hat(3) - param.TP(2)))^2);
    condition = (radius < param.tTol) &...
                (abs(x_hat(2)) < param.rTol) &...
                (abs(x_hat(4)) < param.rTol)&...
                (abs(x_hat(6)) < param.rTol)&...
                (abs(x_hat(8)) < param.rTol);
%     condition = 0; % Leave the condition to be handled by target gen
    if ~condition    
        %% MPC Controller
        opt = mpcqpsolverOptions;
        opt.IntegrityChecks = false;%% for code generation
        opt.FeasibilityTol = 1e-3;
        opt.DataType = 'double';
        %% your code starts here
        % Cholksey and inverse already computed and stored in H
        w = x_hat(1:8) - r(1:8);
        if param.soft == 0
            f = w'*param.G'; % Hard
        else            
            f = [w'*param.G', param.gs']; % Soft
        end
        b = -(param.bb + param.J*x_hat(1:8) + param.L*r(1:8));
        [ubar, ~, ~, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), false(size(param.bb)), opt);
        %% your remaining code here
        if param.mod == 1
%             ubar = param.M1*param.M2*x_hat(1:8) + param.M1*ubar;
            v = param.K_lqr*x_hat(1:8) + ubar(1:2);
        else
            v = ubar(1:2);
        end
        u = v;
    end
end % End of myMPController
