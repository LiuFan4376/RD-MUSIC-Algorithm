% 2D-MUSIC算法

clc; 
clear ; 
close ;

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
K = 300;   %快拍数 
target=2;% 目标数目
Ru=c/(2*delta_f);
theta=(-90:1:90)*pi/180; %测量角度向量
R=linspace(0,Ru,Ru); %测量距离向量

%% ---------------------产生噪声信号
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %产生噪声信号 



%% ------------------------------------------- 产生第一个目标信号
SNR = 10;             
R1 = 1e5;          %delta_r为0
theta1 = 60/180*pi;                          %角度
a1=steer_vector(f0,Delta_f,Dt,Dr,theta1,R1);  %目标一的导向矢量
sig = a1 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)是散射系数
sig1 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 控制信噪比




%% ------------------------------------------- 产生第二个目标信号
SNR = 10;             
R2 = 0.75e5;          %delta_r为0
theta2 = 30/180*pi;                          %角度
a2=steer_vector(f0,Delta_f,Dt,Dr,theta2,R2);  %目标一的导向矢量
sig = a2 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)是散射系数
sig2 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 控制信噪比




%% --------------------------------------------接收信号处理
data = sig1 +sig2+ noise;        %接收数据
% re = res_r*r_num;          % 先验估计距离
% g = kron(ones(N,1),exp(j*2*pi*2*delta_f/c*re*(0:M-1)'));        %补偿矢量
% data_comp = zeros(size(data));         
% for k = 1 : K
%     data_comp(:,k) = data(:,k).*g;             % 补偿后的接收数据
% end
% data_FDA_MIMO =data_comp ;   

Rxx=data*data'/K;

%% 特征值分解
[EV,D]=eig(Rxx);                   %特征值分解
EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
[EVA,I]=sort(EVA);                 %将特征值排序 从小到大
EV=fliplr(EV(:,I));                % 对应特征矢量排序
En=EV(:,target+1:M*N);                  %噪声子空间

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
% [x,y]=find(P==max(max(P)));
% esti_angle=y-91; %算法估计角度
% esti_distance=Ru/1000*(x-1); %算法估计距离
% disp(['target 1: angle =',esti_angle,' distance =',esti_distance]);


figure(1); 
imagesc(theta*180/pi,R,10*log(abs(P)/max(max(abs(P))))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('2D MUSIC算法估计结果');



