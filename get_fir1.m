function [d,d1]=get_fir1(np,f)
% code by He Liu (aresmiki@163.com)
d=[];
d1=[];
    for i=1:length(f)-1
        if i==1
            temp=fir1(np,f(i+1),'low');
        elseif i==length(f)-1
            temp=fir1(np-1,f(i),'high');
            temp=[temp,0]; %≤π¡„
        else
             temp=fir1(np,f(i:i+1),'bandpass');
        end
        d=[d,temp(:)./norm(temp(:))];
        d1=[d1,temp(:)];
    end
    
end
