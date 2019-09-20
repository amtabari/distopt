function x = xModifiedGA (x)
%%%
for i=1:numel(x)
    x(i)=round(x(i));
end

%%% For Genetic Algorithm
for k=1:size(x,1)
    if x(k,1) <= x(k,2)
        x(k,2)=x(k,1)-1;
%         x(k,2)=2;
    end
end

%%% For Genetic Algorithm in 3 Dimensions
for k=1:size(x,1)
    if x(k,3) >= 50
        x(k,3)=2;
    else
        x(k,3)=1;
    end
end
