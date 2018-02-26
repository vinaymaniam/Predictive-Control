function [u, status, iA] = genMPController_Soft(H,G,gs,F,bb,J,L,x,xTarget,m,iA)

% [u,status,iA] = genMPController(H,[(x-xTarget)'*G', gs'],F,bb,J,L,x,xTarget,m,iA);

opt = mpcqpsolverOptions;
opt.IntegrityChecks = false;%% for code generation
opt.FeasibilityTol = 1e-3;
opt.DataType = 'double';
%% This version doesn't make much sense but it works
%% your code starts here
% Cholksey and inverse already computed and stored in H
Linv = H;
w = x - xTarget;
f = [w'*G', gs'];

b = -(bb + J*x + L*xTarget);
[v, status, iA, ~] = mpcqpsolver(Linv, f', -F, b, [], zeros(0,1), iA, opt);
%% your remaining code here
u = v(1:m);


end