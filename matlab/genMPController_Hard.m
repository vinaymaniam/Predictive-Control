function [u, status, iA] = genMPController_Hard(H,G,F,bb,J,L,x,xTarget,m,iA)

[u,status,iA] = genMPController(H,G,F,bb,J,L,x,xTarget,m,iA);

end