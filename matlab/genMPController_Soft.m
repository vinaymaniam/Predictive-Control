function [u, status, iA] = genMPController_Soft(H,G,gs,F,bb,J,L,x,xTarget,m,iA)

[u,status,iA] = genMPController(H,[G; gs],F,b,J,L,x,xTarget,m,iA);

end