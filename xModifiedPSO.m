function x = xModifiedPSO (x,lb,up)
%%%
for i=1:numel(x)
    x(i)=round(x(i));
end


%%% For PSO in 2 Dimensions
for i=1:size(x,2)
    for j=1:2
        if x(j,i)<lb || x(j,i)>up
            x(j,i)=round(lb+rand*(up-lb));
        end
    end
end
for k=2:3:numel(x)
    if x(k) >= x(k-1)
        x(k)=x(k-1)-1;
%         x(k)=2;
    end
end


%%% For PSO in 3 Dimensions
for i=1:size(x,2)
    if x(3,i)~=1 || x(3,i)~=2
        x(3,i)=round(1+rand);
    end
end
