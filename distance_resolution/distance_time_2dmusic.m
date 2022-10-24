clc; 
clear ; 
close ;
% % 时间已过 131.723541 秒。
% % 时间已过 33.498801 秒。
% % 时间已过 14.670724 秒。
% % 时间已过 8.249863 秒。
% % 时间已过 5.492221 秒。
% % 时间已过 4.098453 秒。
% % 时间已过 2.711773 秒。
% % 时间已过 2.066746 秒。
% % 时间已过 1.693376 秒。
% % 时间已过 1.367774 秒
%% ---------------------------参数设置
j = sqrt(-1);
c = 3e8;
M= 8;             % 发射阵元数
N = 6;             % 接收阵元数
f0 = 2e10;          % 参考频率
lamda0 = c/f0;
dt = lamda0/2;      % 发射阵元间距
dr = lamda0/2;      % 接收阵元间距
Dt=(0:M-1)*dt;      % 发射阵列阵元间距设置
Dr=(0:N-1)*dr;      % 接收阵列阵元间距设置 
delta_f= 1000;   % 频率步进量
Delta_f=(0:M-1)*delta_f; %发射阵列频偏设置
K = 500;   %快拍数 
target_num=1;% 目标数目
target_angle=[0];  %目标角度
target_distance=[0.5e3];   % 目标距离

distance_resolution=1:1:10; %距离维分辨率
theta=(-90:1:90)*pi/180;


%%  噪声设置
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; 


%% 目标信号设置
snr=10;
data=noise;  
for i=1:target_num          
    a=steer_vector(f0,Delta_f,Dt,Dr,target_angle(i)*pi/180,target_distance(i));  %导向矢量
    sig = a * sqrt(10^(snr/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)
    sig= sqrt(10^(snr/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 
    data=data+sig;
end


%% --------------------------------------------接收信号处理

for l=distance_resolution
 tic;
  R=1:l:2e4;
Rxx=data*data'/K;

%% 特征值分解
[EV,D]=eig(Rxx);                   %特征值分解
EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
[EVA,I]=sort(EVA);                 %将特征值排序 从小到大
EV=fliplr(EV(:,I));                % 对应特征矢量排序
En=EV(:,target_num+1:M*N);                  %噪声子空间



theta=(-90:l:90)*pi/180; 
    
%% 谱峰搜索
P= zeros(length(theta),length(R)); %波束方向图

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(n),R(m)); %导向矢量
         J=a'*En*En'*a;
        P(n,m) =1/J;
    end
 end
 
P=P';
[x,y]=find(P==max(max(P)));
% esti_angle=y-91; %算法估计角度
% esti_distance=Ru/1000*(x-1); %算法估计距离
% disp(['target 1: angle =',esti_angle,' distance =',esti_distance]);

 toc;
end
% 
% figure(1); 
% imagesc(theta*180/pi,R,10*log(abs(P)/max(max(abs(P))))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('2D MUSIC算法估计结果');
