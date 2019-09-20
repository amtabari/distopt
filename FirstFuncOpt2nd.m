function [fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx,e]=FirstFuncOpt2nd(formul,ax,spp,varb,np,e)
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------

% creation of " mx " matrix , in a random manner
for j=1:spp
    for i=1:varb
        mxrand(i,j,1)=round(1+rand*(np(i)-1));
    end
end
mxrand(:,:,1)=mxrand(:,:,1);

% substitution of created values by " mx " matrix , in " ax " matrix
for j=1:spp
    for i=1:varb
        xrand(i,j)=ax(i,mxrand(i,j));
    end
end
xrand(:,:);

% evaluation of " x " matrix by f @handle , namely " func.m "
for j=1:spp
    fxrand(1,j)=feval(@funcformul,formul,xrand(:,j));
end
fxrand;
e=e+spp;
%----------------------------   Quick  Selection     ---------------------

[minfxsel,pos]=max(fxrand(:));                     % "max" is intentionaly
bestminmxsel(:,1)=mxrand(:,pos);
bestminxsel(:,1)=xrand(:,pos);
bestminfxsel=fxrand(pos);

%-------------------------------   Quick  Total     ----------------------
mxtotal=mxrand;
xtotal =xrand ;
fxtotal=fxrand;
%-------------------------------------------------------------------------
bestfx=bestminfxsel(1);
bestmx=bestminmxsel(:,1);
bestx=bestminxsel(:,1);

