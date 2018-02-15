function [Gamma,Phi] = genPrediction(A,B,N)
% GENPREDICTION  [Gamma,Phi] = genPrediction(A,B,N). 
% A and B are discrete-time state space matrices for x[k+1]=Ax[k]+Bu[k]
% N is the horizon length. 
% Your code is suppose to work for any linear system, not just the gantry crane. 

% Write your code here
n = size(A,1);
E = [A; zeros((N-1)*n,size(A,2))];
A_bar = eye(N*n) - [zeros(n,N*n); [kron(eye(N-1),A), zeros((N-1)*n, n)]];
Phi = A_bar\E;
B_bar = kron(eye(N),B);
Gamma = A_bar\B_bar;

end