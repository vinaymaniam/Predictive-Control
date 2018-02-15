function [DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N)
I = eye(N);

DD = kron(I, Dt);
EE = kron(I, Et);
bb = kron(ones(N,1), bt);
end