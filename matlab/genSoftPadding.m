function [Hs,gs,Fs,bs,Js,Ls] = genSoftPadding(H,F,bb,J,L,S,rho,m)
% S = weight for quadratic cost of constraint violations
% rho = scalar weight for 1-norm of constraint violations
% m = number of inputs
% your code goes here

N = size(H,1)/m;
% H = H(1:m, 1:m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IkronS = 2*kron(eye(N),S);
O1 = zeros(size(H,1), size(IkronS,2));
O2 = zeros(size(IkronS,1), size(H,2));
Hs = [H, O1;
      O2, 2*kron(eye(N),S)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gs = rho*ones(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nConstr = size(S,1);
It = [-eye(nConstr); -eye(nConstr); zeros(m,nConstr); zeros(m,nConstr)];
Ibar = kron(eye(N), It);

Fs = [F, Ibar;
      zeros(size(Ibar,2), size(F,2)), -eye(size(Ibar,2))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs = [bb; zeros(size(Ibar,2), size(bb,2))];
Js = [J; zeros(size(Ibar,2), size(J,2))];
Ls = [L; zeros(size(Ibar,2), size(L,2))]; % L = 0 since ue = 0

end