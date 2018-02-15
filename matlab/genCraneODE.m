function [A,B,C,D] = genCraneODE(m,M,MR,r,g,Tx,Ty,Vm,Ts)
% Inputs:
% m = Pendulum mass (kg)
% M = Cart mass (kg)
% MR = Rail mass (kg)
% r = String length (m)
% g = gravitational accelaration (m*s^-2)
% Tx = Damping coefficient in X direction (N*s*m^-1)
% Ty = Damping coefficient in Y direction (N*s*m^-1)
% Vm = Input multiplier (scalar)
% Ts = Sample time of the discrete-time system (s)
% Outputs:
% A,B,C,D = State Space matrices of a discrete-time or continuous-time state space model

% The motors in use on the gantry crane are identical and therefore Vx=Vy.
Vx=Vm;
Vy=Vm;

% replace A,B,C,D with the correct values
% A=eye(2);
% B=eye(2);
C=eye(8);
D=zeros(8,2);
MMR = M+MR;
AandB = [0,1,0,0,0,0,0,0,0,0;
     0,-1*(Tx/MMR),0,0,g*m/MMR,0,0,0,Vx/MMR,0;
     0,0,0,1,0,0,0,0,0,0;
     0,0,0,-Ty/M,0,0,g*m/M,0,0,Vy/M;
     0,0,0,0,0,1,0,0,0,0;
     0,Tx/(r*MMR),0,0,-1*g*(MMR+m)/(r*MMR),0,0,0,-Vx/(r*MMR),0;
     0,0,0,0,0,0,0,1,0,0;
     0,0,0,Ty/(M*r),0,0,-g*(M+m)/(M*r),0,0,-Vy/(M*r)];
 A = AandB(1:8, 1:8);
 B = AandB(1:8, 9:10);
% if Ts>0 then sample the model with a zero-order hold (piecewise constant) input, otherwise return a continuous-time model
if Ts>0
%     Return discrete-time SS model matrices
    ct = ss(A,B,C,D);
    dt = c2d(ct, Ts);
    A = dt.A;
    B = dt.B;
    C = dt.C;
    D = dt.D;
end

end