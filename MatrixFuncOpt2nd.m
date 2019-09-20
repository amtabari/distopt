function[bestfx,bestmx,bestx,e]=MatrixFuncOpt2nd(formul,varb,nc,np,randomno,NOP,N,rp,...
    lowerbound,upperbound,e,ax,fxtotal,mxtotal,xtotal,bestmx,bestx,bestfx)

    w=1;
    tollfirst=-10;
    while tollfirst < 0
% %                 for loopfirst=1:4
        w=w+1;
         %%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%
 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%     Selection     %%%%%%%%%%%%%%%%%%%
 
        [fxtotalsort,pos]=sort(fxtotal(1,:));
        bestminfxsel(w)=fxtotalsort(1);
        bestminmxsel(:,w)=mxtotal(:,pos(1));
        bestminxsel(:,w) =xtotal(:,pos(1));

% 
        minmxsel=mxtotal(:,pos(1:end));
        minxsel = xtotal(:,pos(1:end));
        minfxsel=fxtotalsort(1:end);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Randomizer " selection "
% 
        for j=1:varb
            for k=1:randomno
                for i=1:varb
                    mxrandsel(i,k,j)=round(1+rand*(np(i)-1));
                end
                mxrandsel(j,:,j)=bestminmxsel(j,w);
            end
        end
        mxrandsel;

        % substitution of created values by " mx " matrix , in " ax " matrix
        for z=1:varb
            for j=1:randomno%length(mxrandsel)
                for i=1:varb
                    xrandsel(i,j,z)=ax(i,mxrandsel(i,j,z));
                end
            end
            
            % evaluation of " x " matrix by f @handle , namely " func.m "
            for j=1:randomno%length(mxrandsel)
                fxrandsel(1,j,z)=feval(@funcformul,formul,xrandsel(:,j,z));
            end
            e=e+randomno;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fxrandselsort,pos]=sort(fxrandsel(1,:));
        bestminfxrandsel=fxrandselsort(1);
        minfxrandsel=fxrandselsort(1:randomno);

        minmxrandsel=mxrandsel(:,pos(1:randomno));
        bestminmxrandsel=minmxrandsel(:,1);

        minxrandsel=xrandsel(:,pos(1:randomno));
        bestminxrandsel=minxrandsel(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[minfxsel,minfxrandsel];
        [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[minmxsel,minmxrandsel];
        mxtotal=mxtotal(:,pos(1:end));
        bestmxtotal=mxtotal(:,1);

        xtotal=[minxsel,minxrandsel];
        xtotal=xtotal(:,pos(1:end));
        bestxtotal=xtotal(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Rain Algorithm    %%%%%%%%%%%%%%

        for p=1:varb

            x=linspace(lowerbound(p),upperbound(p),NOP);
            y=zeros;
            for j=1:NOP
                xforward=bestxtotal(:,end);
                xforward(p)=x(j);
                
                y(j)=feval(@funcformul,formul,xforward);
            end
            e=e+NOP;

            [yvall(2,1),pos]=min(y(:));
            vall(2,1)=x(pos);
            v=0;
            for i=pos:-1:1
                if y(i) > yvall(2,1)
                    yvall(1,1)=y(i);
                    vall(1,1)=x(i);
                    v=v+0.5;
                    break
                end
            end
            for i=pos:NOP
                if y(i) > yvall(2,1)
                    yvall(3,1)=y(i);
                    vall(3,1)=x(i);
                    v=v+0.5;
                    break
                end
            end

            if v ~= 1
                v=0;
            end
            for i=2:(NOP-1)
                if y(i)<y(i+1) && y(i)<y(i-1)
                    v=v+1;

                    vall(1,v)=x(i-1);
                    yvall(1,v)=y(i-1);
                    vall(2,v)=x(i);
                    yvall(2,v)=y(i);
                    vall(3,v)=x(i+1);
                    yvall(3,v)=y(i+1);
                end
            end
            if v == 0
                [Fxminimum(p),pos]=min([ y(1) , y(end) ]);
                if pos==1
                    Xminimum(p)=x(1);
                else
                    Xminimum(p)=x(NOP);
                end

            end
            if v >= 1

                Fxmin=zeros;
                for valley=1:v

                    a=vall(1,valley);
                    b=vall(3,valley);

                    x=linspace(a,b,N);
                    y=zeros;
                    for j=1:N
                        xforward=bestminxsel(:,end);
                        xforward(p)=x(j);
                        
                        y(j)=feval(@funcformul,formul,xforward);
                    end
                    e=e+ N ;

                    xmm=zeros;
                    xmm(1)=vall(2,valley);
                    k=1;
                    for j=1:(N-3)
                        for i=(j+1):(N-2)
                            m1=(y(j+1)-y(j))/(x(j+1)-x(j));
                            m2=(y(i+2)-y(i+1))/(x(i+2)-x(i+1));
                            intersection=(-y(i+1)+m2*x(i+1)+y(j)-m1*x(j))/(m2-m1) ;
                            if intersection > a && intersection < b
                                k=k+1;
                                xmm(k)=intersection;
                            end
                            for r=1:rp
                                k=k+1;
                                xmm(k)=vall(1,valley)+rand*(vall(3,valley)-vall(1,valley));
                            end
                        end
                    end
                    fxmm=zeros;
                    for j=1:numel(xmm)
                        xmmforward=bestminxsel(:,end);
                        xmmforward(p)=xmm(j);
                        
                        fxmm(j)=feval(@funcformul,formul,xmmforward);
                    end
                    e=e+ numel(xmm) ;
                    [Fxmin(valley),pos]=min(fxmm(:));
                    Xmin(valley)=xmm(pos);

                end
                [Fxminimum(p),pos]=min(Fxmin(:));
                Xminimum(p)=Xmin(pos);
                Fxminimum(p);

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%    Finding nearest array
            %%%%%%%%%%%%%%%%%%%%%%%%%%%      in  " ax Matrix "
            dist=zeros;

            for j=1:np
                dist(j)=abs(ax(p,j)-Xminimum(p)) ;   % ax should be changed
            end

            minxrain(:,p)=bestxtotal(:,end);
            [mindist,pos]=min(dist(:));
            minxrain(p,p)=ax(p,pos);                 % ax should be changed
            
            minfxrain(p)=feval(@funcformul,formul,minxrain(:,p));
            e=e+ 1 ;

            minmxrain(:,p)=bestmxtotal(:,end);
            minmxrain(p,p)=pos;
        end

        [minfxrainsort,pos]=sort(minfxrain(1,:));
        bestminfxrain=minfxrainsort(1);
        minfxrain=minfxrainsort(1:varb);

        minmxrain=minmxrain(:,pos(1:varb));
        bestminmxrain=minmxrain(:,1);

        minxrain =minxrain(:,pos(1:varb));
        bestminxrain=minxrain(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[fxtotal,minfxrain];
        [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[mxtotal,minmxrain];
        mxtotal=mxtotal(:,pos(1:end));
        bestmxtotal=mxtotal(:,1);

        xtotal=[xtotal,minxrain];
        xtotal=xtotal(:,pos(1:end));
        bestxtotal=xtotal(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Randomizer " Rain "
        randomno;
        for j=1:varb
            for k=1:randomno
                for i=1:varb
                    mxrandrain(i,k,j)=round(1+rand*(np(i)-1));
                end
                mxrandrain(j,:,j)=bestmxtotal(j);
            end
        end
        mxrandrain;

        % substitution of created values by " mx " matrix , in " ax " matrix
        for z=1:varb
            for j=1:randomno%length(mxrandrain)
                for i=1:varb
                    xrandrain(i,j,z)=ax(i,mxrandrain(i,j,z));
                end
            end
           
            % evaluation of " x " matrix by f @handle , namely " func.m "
            for j=1:randomno%length(mxrandrain)
                fxrandrain(1,j,z)=feval(@funcformul,formul,xrandrain(:,j,z));
            end
            e=e+ randomno ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fxrandrainsort,pos]=sort(fxrandrain(1,:));
        bestminfxrandrain=fxrandrainsort(1);
        minfxrandrain=fxrandrainsort(1:randomno);

        minmxrandrain=mxrandrain(:,pos(1:randomno));
        bestminmxrandrain=minmxrandrain(:,1);

        minxrandrain=xrandrain(:,pos(1:randomno));
        bestminxrandrain=minxrandrain(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[fxtotal,minfxrandrain];
        [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[mxtotal,minmxrandrain];
        mxtotal=mxtotal(:,pos(1:end));
        bestmxtotal=mxtotal(:,1);

        xtotal=[xtotal,minxrandrain];
        xtotal=xtotal(:,pos(1:end));
        bestxtotal=xtotal(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Combination Old
        %     for k=1:nc
        %         for j=1:2:(fix(length(mxtotal)/2))*2
        %             for i=1:k
        %                 mxcomb(i,j,k)=mxtotal(i,j+1);
        %                 mxcomb(i,j+1,k)=mxtotal(i,j);
        %             end
        %             for i=(k+1):varb
        %                 mxcomb(i,j,k)=mxtotal(i,j);
        %                 mxcomb(i,j+1,k)=mxtotal(i,j+1);
        %             end
        %         end
        %     end
        %     mxcomb;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Combination New
        for k=1:nc
            for j=1:length(mxtotal)
                for i=1:(k+1)
                    mxcomb(i,j,k)=bestmxtotal(i,1);
                end
                for i=(k+2):varb
                    mxcomb(i,j,k)=mxtotal(i,j);
                end
            end
        end
        mxcomb;
        % substitution of created values by " mx " matrix , in " ax " matrix
        for z=1:(k)
            for j=1:length(mxcomb)
                for i=1:varb
                    xcomb(i,j,z)=ax(i,mxcomb(i,j,z));
                end
            end
        
            % evaluation of " x " matrix by f @handle , namely " func.m "
            for j=1:length(mxcomb)
                fxcomb(1,j,z)=feval(@funcformul,formul,xcomb(:,j,z));
            end
            e=e+ length(mxcomb) ;
        end
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fxcombsort,pos]=sort(fxcomb(1,:));
        bestminfxcomb=fxcombsort(1);
        minfxcomb=fxcombsort(1:end);

        minmxcomb=mxcomb(:,pos(1:end));
        bestminmxcomb=minmxcomb(:,1);

        minxcomb=xcomb(:,pos(1:end));
        bestminxcomb=minxcomb(:,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[fxtotal,minfxcomb];
        [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[mxtotal,minmxcomb];
        mxtotal=mxtotal(:,pos(1:end));
        bestmxtotal=mxtotal(:,1);

        xtotal=[xtotal,minxcomb];
        xtotal=xtotal(:,pos(1:end));
        bestxtotal=xtotal(:,1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Randomizer " first array "
        randomno;
        for j=1:varb
            for k=1:randomno
                for i=1:varb
                    mxrandcomb(i,k,j)=round(1+rand*(np(i)-1));
                end
                mxrandcomb(j,:,j)=bestmxtotal(j);
            end
        end
        mxrandcomb;

        % substitution of created values by " mx " matrix , in " ax " matrix
        for z=1:varb
            for j=1:randomno%length(mxrandcomb)
                for i=1:varb
                    xrandcomb(i,j,z)=ax(i,mxrandcomb(i,j,z));
                end
            end
         
            % evaluation of " x " matrix by f @handle , namely " func.m "
            for j=1:randomno%length(mxrandcomb)
                fxrandcomb(1,j,z)=feval(@funcformul,formul,xrandcomb(:,j,z));
            end
            e=e+ randomno ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fxrandcombsort,pos]=sort(fxrandcomb(1,:));
        bestminfxrandcomb=fxrandcombsort(1);
        minfxrandcomb=fxrandcombsort(1:randomno);

        minmxrandcomb=mxrandcomb(:,pos(1:randomno));
        bestminmxrandcomb=minmxrandcomb(:,1);

        minxrandcomb=xrandcomb(:,pos(1:randomno));
        bestminxrandcomb=minxrandcomb(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[fxtotal,minfxrandcomb];
        [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[mxtotal,minmxrandcomb];
        mxtotal=mxtotal(:,pos(1:end));
        bestmxtotal=mxtotal(:,1);

        xtotal=[xtotal,minxrandcomb];
        xtotal=xtotal(:,pos(1:end));
        bestxtotal=xtotal(:,1);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Randomizer " second array "
        randomno;
        for j=1:(varb-1)
            for k=1:randomno
                for i=1:varb
                    mxrandsecond(i,k,j)=round(1+rand*(np(i)-1));
                end
            end
            for m=j:j+1
                mxrandsecond(m,:,j)=bestmxtotal(m);
            end
        end
        mxrandsecond;

        % substitution of created values by " mx " matrix , in " ax " matrix
        for z=1:(varb-1)
            for j=1:randomno%length(mxrandsecond)
                for i=1:varb
                    xrandsecond(i,j,z)=ax(i,mxrandsecond(i,j,z));
                end
            end
           
            % evaluation of " x " matrix by f @handle , namely " func.m "
            for j=1:randomno%length(mxrandsecond)
                fxrandsecond(1,j,z)=feval(@funcformul,formul,xrandsecond(:,j,z));
            end
            e=e+ randomno ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fxrandsecondsort,pos]=sort(fxrandsecond(1,:));
        bestminfxrandsecond=fxrandsecondsort(1);
        minfxrandsecond=fxrandsecondsort(1:end);

        minmxrandsecond=mxrandsecond(:,pos(1:end));
        bestminmxrandsecond=minmxrandsecond(:,1);

        minxrandsecond=xrandsecond(:,pos(1:end));
        bestminxrandsecond=minxrandsecond(:,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Total Values
        fxtotal=[bestminfxsel(end),bestminfxrandsel,bestminfxrain,...
            bestminfxrandrain,bestminfxcomb,bestminfxrandcomb,bestminfxrandsecond];
        %     [fxtotal,pos]=sort(fxtotal(1,:));

        mxtotal=[bestminmxsel(:,end),bestminmxrandsel,bestminmxrain,...
            bestminmxrandrain,bestminmxcomb,bestminmxrandcomb,bestminmxrandsecond];
        %     mxtotal=mxtotal(:,pos(1:end));
        xtotal=[bestminxsel(:,end),bestminxrandsel,bestminxrain,...
            bestminxrandrain,bestminxcomb,bestminxrandcomb,bestminxrandsecond];
        %     xtotal=xtotal(:,pos(1:end));


        %%%%%%%%%%%%%%%%%%%%%    Tolerance First
        tollfirst=bestminfxsel(w-1)-bestminfxsel(w);
    end
    bestfx=bestminfxsel(end);
    bestmx(:)=bestminmxsel(:,end);
    bestx(:)=bestminxsel(:,end);

%   w
% end


