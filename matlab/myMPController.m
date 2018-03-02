function u = myMPController(r, x_hat, param)
    %% Do not delete this line
    % Create the output array of the appropriate size
    u = zeros(2,1);
    %% MPC Controller
    opt = mpcqpsolverOptions;
    opt.IntegrityChecks = false;%% for code generation
    opt.FeasibilityTol = 1e-3;
    opt.DataType = 'double';
    %% This version doesn't make much sense but it works
    %% your code starts here
    % Cholksey and inverse already computed and stored in H
    w = x_hat - r;
    f = w'*param.G';

    b = -(param.bb + param.J*x_hat + param.L*r);
    [v, ~, ~, ~] = mpcqpsolver(param.H, f', -param.F, b, [], zeros(0,1), false(size(param.bb)), opt);
    %% your remaining code here
    u = v(1:2);

end % End of myMPController
