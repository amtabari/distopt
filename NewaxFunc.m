function[axnew]=NewaxFunc(t,varb,bestmx,secondp,ax,np)

t=1;
    for i=1:varb
        if bestmx(i,t) == 1          % if  " bestmx "  was  " 1 "
            axnew(i,1:(2+2*secondp)) = linspace(ax(i,1),ax(i,2),(2+2*secondp));
            axnew(i,(3+2*secondp):(np(i)+2*secondp)) = ax(i,3:np(i));
        elseif bestmx(i,t) == np(i)  % if  " bestmx "  was  " last "
            axnew(i,(np(i)-1):(np(i)+2*secondp)) = linspace(ax(i,(np(i)-1)),ax(i,np(i)),(2+2*secondp));
            axnew(i,1:(np(i)-2)) = ax(i,1:(np(i)-2));
        else                                 % if  " bestmx "  was  Normal
            intervalleft=linspace(ax(i,bestmx(i,t)-1),ax(i,bestmx(i,t)),(2+secondp));
            axnew(i,(bestmx(i,t)-1):((bestmx(i,t)-1)+secondp)) = intervalleft(1:1+secondp);

            intervalright=linspace(ax(i,bestmx(i,t)),ax(i,bestmx(i,t)+1),(2+secondp));
            axnew(i, (((bestmx(i,t)-1)+secondp)+1) : (((bestmx(i,t)-1)+secondp)+numel(intervalright)) ) = intervalright(1:end);

            axnew(i,1:bestmx(i,t)-2) = ax(i,1:bestmx(i,t)-2);
            axnew(i,(bestmx(i,t)+1)+(2*secondp+1):(np(i)+2*secondp)   )=ax(i,(bestmx(i,t)+2):np(i));
        end
    end
    