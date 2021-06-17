function [optWQ,recQ,obj,optWQ2,recQ2,obj2,rec,rec2]=min_lplq_optimation(x,np,nlevel,p,q)
    %  code by Liu He (aremiki@163.com), January 2020
    %  Used in my PhD research at the University of SouthWest Jiaotong University.
    %   
    %
    %    Algorithm Reference:
    %    Optimized minimum generalized Lp/Lq deconvolution for recovering repetitive 
    %    impacts from a vibration mixture, Measurement, vol. 168，
    %     doi: 10.1016/j.measurement.2020.108329 
    % optWQ,recQ,obj 优化后的结果
    % optWQ2,recQ2,obj2 直接滤波的结果
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if nargin <= 3
    p=1;q=2;
end
%%  具有初始化能力的最优追踪
    %% ?初始化滤波系数 
        %% 构造Hankel 矩阵
    N=length(x);
    A=[];
    for iI=1:floor(N-np+1)
        A(iI,:)=x(iI:(iI+np-1));
    end
    Level_w = 1:nlevel;	Level_w = [Level_w;Level_w+log2(3)-1];	Level_w = Level_w(:); Level_w = [0 Level_w(1:2*nlevel-1)'];
    
    %% 塔式
    for i=1:length(Level_w)
        if i==1
            optSWM(:,i)=zeros(np,1);
            optSWM(2,i)=1;
            optSWM1=optSWM;
        else
            if(mod(i,2)==0)
                qi=2^floor(i/2);
            else
                qi=3*2^(floor(i/2)-1);
            end
             f=0:1/qi:1;
             [opt_temp,opt_temp1]=get_fir1(np-1,f);
             optSWM=[optSWM,opt_temp];
             optSWM1=[optSWM1,opt_temp1];
        end
    end
    for i=1:size(optSWM,2)
        [optW(:,i),rec(:,i),obj(i)]=min_L1_L2(A,0,optSWM(:,i),p,q);
       rec2(:,i)=A*optSWM(:,i);
       obj2(i)=sign(log(q/p)).*(norm(abs(rec2(:,i)),p)/norm(abs(rec2(:,i)),q)).^p;
    end
    obj=obj./sqrt(N-np+1);
    obj2=obj2./sqrt(N-np+1);
    [val,idx]=min(obj);
    optWQ=optW(:,idx);
    recQ=rec(:,idx);
    [val2,idx2]=min(obj2);
    optWQ2=optSWM(:,idx2);
    recQ2=rec2(:,idx2);
    
    %% 绘图
Fs=10000;
freq_w = Fs*((0:3*2^nlevel-1)/(3*2^(nlevel+1)) + 1/(3*2^(2+nlevel)));
Kwav=zeros(2*nlevel,length(freq_w));
Kwav2=Kwav;
iter_idx=1;
 for i=1:length(Level_w)
        if i==1
           Kwav(i,:)=obj(iter_idx);
           Kwav2(i,:)=obj2(iter_idx);
           iter_idx=iter_idx+1;
        else
            if(mod(i,2)==0)
                qi=2^floor(i/2);
            else
                qi=3*2^(floor(i/2)-1);
            end
            for j=1:qi
                Kwav(i,(j-1)*(length(freq_w)/qi)+1:j*(length(freq_w)/qi))=obj(iter_idx);
                Kwav2(i,(j-1)*(length(freq_w)/qi)+1:j*(length(freq_w)/qi))=obj2(iter_idx);
                iter_idx=iter_idx+1;
            end
        end
 end
Kwav=-Kwav;
Kwav2=-Kwav2;
[I,J,M] = max_IJ(Kwav);
fi = (J-1)/3/2^(nlevel+1);   fi = fi + 2^(-2-Level_w(I));
figure
imagesc(freq_w,1:2*nlevel,Kwav),colorbar,
colormap(jet)
colorbar
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,300,220]);
xlim([0,Fs/2])
xlabel('Frequency (Hz)','fontsize',12),set(gca,'ytick',1:2*nlevel,'yticklabel',round(Level_w*10)/10),ylabel('Level k','fontsize',12)
% title(['fb-kurt.2 - K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'],'fontsize',12)
title(['K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'],'fontsize',12)
%%
[I,J,M] = max_IJ(Kwav2);
fi = (J-1)/3/2^(nlevel+1);   fi = fi + 2^(-2-Level_w(I));
figure
imagesc(freq_w,1:2*nlevel,Kwav2),colorbar,
colormap(jet)
colorbar
set(gca,'linewidth',1);
set(gca,'FontSize',12);
set(gcf,'position',[200,300,300,220]);
xlim([0,Fs/2])
xlabel('Frequency (Hz)','fontsize',12),set(gca,'ytick',1:2*nlevel,'yticklabel',round(Level_w*10)/10),ylabel('Level k','fontsize',12)
% title(['fb-kurt.2 - K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'],'fontsize',12)
title(['K_{max}=',num2str(round(10*M)/10),' @ Level ',num2str(fix(10*Level_w(I))/10),', Bw= ',num2str(Fs*2^-(Level_w(I)+1)),'Hz, f_c=',num2str(Fs*fi),'Hz'],'fontsize',12)

end


function [optW,rec,obj]=min_L1_L2(data,optM,optW,p,q)

%%
if optM==0
[optW,obj] = minFunc(@L1_L2obj_grad, optW(:), ...
                   struct('MaxIter', 400,'Display','off'), data,p,q);
%     hold off
else
[optW,obj]  = minFunc(@L1_L2obj_grad, optW(:), ...
                   struct('MaxIter', 400,'numDiff',1,'Display','final'), data,p,q);
end
optW=optW./norm(optW);        
rec= data*optW;
end
%%  目标函数
function obj=L1_L2obj(w,A)
%     w=w./norm(w);
    y=A*w;
    obj=norm(y,1)/norm(y,2);
end
function [obj,grad]=L1_L2obj_grad(w,A,p,q)
    y=A*w;
    C= sqrt(y.^2 + 1e-8);
    %% 可用log(q/p) 来测量到底因该最大化还是最小化 log(2/1) 最小化, log(2/4)最大化 
    obj=sign(log(q/p)).*(norm(C,p)./norm(C,q)).^p;
    aJdR=sign(log(q/p)).*p*(norm(C,p)./norm(C,q)).^(p-1);
       %%    链式法则求导数 p=1,q=2
    aRdC=(norm(C,p).^(1-p).*C.^(p-1).*norm(C,q))./norm(C,q).^2-(norm(C,q).^(1-q).*C.^(q-1).*norm(C,p))./norm(C,q).^2;
    aCdy= y./C;
    grad=A'*(aJdR.*aRdC.*aCdy);  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I,J,M] = max_IJ(X)
% Returns the row and column indices of the maximum in matrix X.

[temp,tempI] = max(X);
[M,J] = max(temp);
I = tempI(J);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
