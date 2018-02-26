function [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N)
% disp(size(DD))
% disp(size(EE))
% disp(size(Gamma))
% disp(size(Phi))
% disp(N)
n = size(DD,2)/N;
m = size(EE,2)/N;
Pt = [eye(n); Phi(1:end-n, :)];
Gt = [zeros(n, N*m); Gamma(1:end-n, :)];


F = DD*Gt + EE;
J = -DD*Pt;
L = -J - DD*kron(ones(N,1),eye(size(Pt,2)));
end