function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% MPC Controller
    opt = mpcqpsolverOptions;
    opt.IntegrityChecks = false;%% for code generation
    opt.FeasibilityTol = 1e-3;
    opt.DataType = 'double';    
    %% Check if we crossed the turning point yet
    if x_hat(1) < param.x_star
        w = x_hat(1:8) - r(1:8);
        f = w'*param.G1';
        b = -(param.bb1 + param.J1*x_hat(1:8) + param.L1*r(1:8));
        [v, ~, ~, ~] = mpcqpsolver(param.H1, f', -param.F1, b, [], zeros(0,1), false(size(param.bb1)), opt);        
    else
        w = x_hat(1:8) - r(1:8);
        f = w'*param.G2';
        b = -(param.bb2 + param.J2*x_hat(1:8) + param.L2*r(1:8));
        [v, ~, ~, ~] = mpcqpsolver(param.H2, f', -param.F2, b, [], zeros(0,1), false(size(param.bb2)), opt);        
    end   
    u = v(1:2);
end % End of myMPController
