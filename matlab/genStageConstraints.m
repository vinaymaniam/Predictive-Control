function [Dt,Et,bt]=genStageConstraints(A,B,D,cl,ch,ul,uh)
% sprintf('(n,m)=(%i,%i)',size(A,2),size(B,2))
DA = D*A;
DB = D*B;
I = eye(size(B,2));
O = zeros(size(B,2),size(A,2));
Dt = [DA; -DA; O; O];
Et = [DB; -DB; I; -I];
bt = [ch; -cl; uh; -ul];
end