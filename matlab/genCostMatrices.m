function [H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N)
%% cost function matrices
% Gamma and Phi are the prediction matrices
% Q is the stage cost weight on the states, i.e. x'Qx
% R is the stage cost weight on the inputs, i.e. u'Ru
% P is the terminal weight on the final state

% Your code goes here
Q_bar = [[kron(eye(N-1),Q), zeros((N-1)*size(Q,1),size(P,2))]; [zeros(size(P,1),(N-1)*size(Q,2)), P]];
R_bar = kron(eye(N),R);
H = 2.*(Gamma'*Q_bar*Gamma + R_bar);
G = 2.*(Phi'*Q_bar*Gamma);
G = G';

end