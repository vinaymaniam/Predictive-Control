function [Hs,gs,Fs,bs,Js,Ls] = genSoftPadding(H,F,bb,J,L,S,rho,m)
% S = weight for quadratic cost of constraint violations
% rho = scalar weight for 1-norm of constraint violations
% m = number of inputs
% your code goes here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hs = [H(1:m, 1:m); S];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs = [eye(m), ?]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q4 = [-ones(?, ?), ones(?, ?), zeros(?, ?), zeros(?, ?)];
Fs = [F, Q4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs = bb;
Js = J;
Ls = L;

end