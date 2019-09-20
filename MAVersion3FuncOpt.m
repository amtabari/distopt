function[bestfxnew,bestmxnew,bestxnew,ax,np,e]=MAVersion3FuncOpt(F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,...
    ax,spp,varb,secondp,np,e,nc,randomno,NOP,N,rp,lowerbound,upperbound)
%-------------------------------------------------------------------------


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%                                FirstFunc
[fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx,e]=feval(@FirstFuncOpt,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,ax,spp,varb,np,e);
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%                                MatrixFunc
[bestfx,bestmx,bestx,e]=feval(@MatrixFuncOpt,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,varb,nc,np,randomno,NOP,N,rp,...
    lowerbound,upperbound,e,ax,fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx);
% Change in bestfx
bestfxnew(1)=bestfx;
bestmxnew(:,1)=bestmx;
bestxnew(:,1)=bestx;
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
t=1;
tollfirst=10;
counter=0;
% while tollfirst > 0  %|| tollfirst == 0
% while tollfirst > 0.001 || counter < 3
while tollfirst > 0 || counter < 3

    t=t+1;
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %                                New ax Func            
    [axnew]=feval(@NewaxFunc,t,varb,bestmx,secondp,ax,np);
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    npnew(1:varb)=np(:)+2*secondp;
    %---------------------------------------------------------------------
    np=npnew;
    ax=axnew;
    %---------------------------------------------------------------------
    abc=0;
    %     tollsecond=-5;
    %     while tollsecond < 0 %&& abc < 3 %|| tollsecond == 0
    for j=1:3
        %     for h=1:7
        abc=abc+1;
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %                                FirstFunc
        [fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx,e]=feval(@FirstFuncOpt,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,ax,spp,varb,np,e);
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %                               MatrixFunc
        [bestfx,bestmx,bestx,e]=feval(@MatrixFuncOpt,F,D,R,TF,bldlupper,bhdhlower,Xcomp,fs_uipanel23,varb,nc,np,randomno,NOP,N,rp,...
            lowerbound,upperbound,e,ax,fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx);
        % Change in bestfx
        bestfxnew(t)=bestfx;
        bestmxnew(:,t)=bestmx;
        bestxnew(:,t)=bestx;
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %--------------------------    Tolerance Second    ---------------
        tollsecond=bestfxnew(t-1)-bestfxnew(t);
        if tollsecond < 0
            bestfxnew(t)=bestfxnew(t-1);
            bestmxnew(:,t)= bestmxnew(:,t-1);
            bestxnew(:,t)=bestxnew(:,t-1);
        end
    end
    abc;
    ax;
    bestmx(:,t)=bestmxnew(:,t);
    %--------------------------    Tolerance First    --------------------
    tollfirst=bestfxnew(t-1)-bestfxnew(t);
    if tollfirst > 0
        counter=0;
    end
    if tollfirst==0
        counter=counter+1;
    end
    fff=bestfxnew(t)
    tollfirst=abs(bestfxnew(t-1)-bestfxnew(t));
end

bestfxnew;
bestmxnew;
bestxnew;

