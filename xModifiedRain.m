function x = xModifiedRain (x)
%%%
for i=1:numel(x)
    x(i)=round(x(i));
end



% %%% For Rain Algorithm in 2 Dimensions
% for k=2:2:numel(x)
%     if x(k) >= x(k-1)
%         x(k)=x(k-1)-1;
% %         x(k)=2;
%     end
% end

%%% For Rain Algorithm in 3 Dimensions
for k=2:3:numel(x)
    if x(k) >= x(k-1)
        x(k)=x(k-1)-1;
%         x(k)=2;
    end
end

