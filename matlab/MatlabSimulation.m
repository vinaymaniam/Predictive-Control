%% This script runs a linear simulation of the plant and tests your controller performance
% create waiting bar
hw=waitbar(0,'Please wait...');
warning('on');
% create statespace object object
sysd=ss(A,B,C,0,Ts);
% starting variables
x=x0;
allU=[];
allOpt=[];
allIter=[];
% initial vector for 'cold start'. see mpcqpsolver
iA = false(size(bb));
for t=0:Ts:T
    waitbar(t/T,hw,'Please wait...');
    tic;
    if(hardConstraints)
        [u,status,iA] = genMPController_Hard(H,G,F,bb,J,L,x(:,end),xTarget,size(B,2),iA);
    else
        [u,status,iA] = genMPController_Soft(H,G,gs,F,bb,J,L,x(:,end),xTarget,size(B,2),iA);
    end
    optTime=toc;    
    if status<0
        close(hw);
        warning('QP solver failed to find a solution!');
        break
    end
    % Simulate
    [yy,tt,xx] =lsim(sysd,[u';0 0],[0 Ts],x(:,end));
    % Store variables
    x=[x xx(end,:)'];
    allU=[allU ;u'];
    allOpt=[allOpt;optTime];
    allIter=[allIter;status];
end
t=0:Ts:t;
x=x(:,1:length(t))';
close(hw);
figure('Name','Optimisation time'); 
bar(t,allOpt);
xlabel('Simulation time [s]')
ylabel('Optimisation time [s]')