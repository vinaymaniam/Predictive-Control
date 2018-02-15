n = size(A,1);

E = [A; zeros((N-1)*n,size(A,2))];
A_bar = eye(N*n) - [zeros(n,N*n); [kron(eye(N-1),A), zeros((N-1)*n)]];

Phi = A_bar\E;

B_bar = kron(eye(N),B);

Gamma = A_bar\B_bar;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

Q_bar = [[kron(eye(N-1),Q), zeros((N-1)*size(Q,1),size(P,2))]; [zeros(size(P,1),(N-1)*size(Q,2)), P]];
R_bar = kron(eye(N),R);
H = 2.*(Gamma'*Q_bar*Gamma + R_bar);
G = 2.*(Phi'*Q_bar*Gamma);
G = G';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

L = H\G;
K = -[eye(m), zeros(m,size(H,1)-m)]*L;