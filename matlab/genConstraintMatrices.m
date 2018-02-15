function [F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N)
% disp(size(DD))
% disp(size(EE))
% disp(size(Gamma))
% disp(size(Phi))
% disp(N)
F = DD*Gamma + EE;
J = -DD*Phi;
L = -J - DD*kron(ones(N,1),eye(size(Phi,2)));
end