function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% MPC Controller
    opt = mpcqpsolverOptions;
    opt.IntegrityChecks = false;%% for code generation
    opt.FeasibilityTol = 1e-3;
    opt.DataType = 'double';    
    persistent iA;
    persistent iA2;
    if isempty(iA)
        iA = false(size(param.bb2));
        iA2 = false(size(param.bb2));
    end    
    radius = sqrt((abs(x_hat(1) - param.TP2(1)))^2+(abs(x_hat(3) - param.TP2(2)))^2);
    condition = (radius < param.tTol^2) &...
                (abs(x_hat(5)) < param.rTol)&...
                (abs(x_hat(7)) < param.rTol);
    if ~condition
        %% Check if we crossed the turning point yet
        if r(1,1) == param.TP1(1) && r(3,1) == param.TP1(2)
            w = x_hat(1:8) - r(1:8);
            f = w'*param.G1';
            b = -(param.bb1 + param.J1*x_hat(1:8) + param.L1*r(1:8));
            [v, ~, ~, ~] = mpcqpsolver(param.H1, f', -param.F1, b, [], zeros(0,1), false(size(param.bb1)), opt);        
            u = v(1:2);
        else
            w = x_hat(1:8) - r(1:8);
            f = w'*param.G2';
            b = -(param.bb2 + param.J2*x_hat(1:8) + param.L2*r(1:8));

            [v, ~, iA2, ~] = mpcqpsolver(param.H2, f', -param.F2, b, [], zeros(0,1), iA2, opt);        
            u = v(1:2);
        end
    end
end % End of myMPController
