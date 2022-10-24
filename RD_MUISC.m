%% RD-MUSIC算法
clc;
close;
clear;

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
target=3;% 目标数目
Ru=c/(2*delta_f);
theta=(-90:0.1:90)*pi/180; %测量角度向量
R=linspace(0,Ru,1000); %测量距离向量

%% ---------------------产生噪声信号
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %产生噪声信号 

%% ------------------------------------------- 产生第一个目标信号
SNR = 10;             
R1 = 10e4;          %delta_r为0
theta1 = -20/180*pi;                          %角度
a1=steer_vector(f0,Delta_f,Dt,Dr,theta1,R1);  %目标一的导向矢量
sig = a1 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)是散射系数
sig1 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 控制信噪比

%% ------------------------------------------- 产生第二个目标信号
SNR = 10;             
R2 = 7.5e4;          %delta_r为0
theta2 = -10/180*pi;                          %角度
a2=steer_vector(f0,Delta_f,Dt,Dr,theta2,R2);  %目标一的导向矢量
sig = a2 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)是散射系数
sig2 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 控制信噪比
SNR = 10;
%% ------------------------------------------- 产生第三个目标信号
SNR = 10;
R3 = 8.5e4;          %delta_r为0
theta3 = 20/180*pi;                          %角度
a3=steer_vector(f0,Delta_f,Dt,Dr,theta3,R3);  %目标一的导向矢量
sig = a3 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)是散射系数
sig3 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 控制信噪比

%% --------------------------------------------接收信号处理
data = sig1 +sig2+sig3+noise;        %接收数据  
 Rxx=data*data'/K;                %协方差矩阵

%% 特征值分解
% [EV,D]=eig(Rxx);                   %特征值分解
% EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
% [EVA,I]=sort(EVA);                 %将特征值排序 从小到大
% EV=fliplr(EV(:,I));                % 对应特征矢量排序
% En=EV(:,target+1:M*N);                  %噪声子空间
G=Rxx(:,1:target);
H=Rxx(:,target+1:M*N);
P=inv(G'*G)*G'*H;
I=-diag(ones(1,M*N-target));
En=[P;I];
En=orth(En);
En = En./repmat(sqrt(sum(En.^2,1)),size(En,1),1);


%% -----------------------------------------------------角度估计
P= zeros(1,length(theta)); %波束方向图
 for n = 1 : length(theta)
         
        d = exp(j*2*pi*f0/c*Dt'*sin(theta(n)));        %  发射阵列角度导向矢量
         b = exp(j*2*pi*f0/c*Dr'*sin(theta(n)));        %  接收阵列导向矢量 
          W=kron(b,diag(d))'*En*En'*kron(b,diag(d));
          W1=W(1,1);
          W2=W(1,2:M);
          W4=W(2:M,2:M);
          J=W1-W2/W4*W2';
        P(n) =1/J;
   
 end
%  P=P';
figure(1);
P=abs(P)/max(abs(P));
plot(theta*180/pi,10*log10(P));
xlabel('角度/o');  ylabel('归一化幅度/dB'); 

axis([-90.1 90.1 -35 0]);
%% 检索角度估计值
estimate_theta=zeros(1,target); %RD-MUSIC算法角度维估计结果
flag=1;
for i=1:length(P)
    if flag==target+1
        break;
    end
    if P(i)<=0.01
        continue;
    end
 
    if i==1&&P(i)>P(i+1)
        estimate_theta(flag)=i;
        flag=flag+1;
    end
    if i==length(P)&&P(i)>P(i-1)
        estimate_theta(flag)=i;
        flag=flag+1;
    end
   if i<length(P)&&i>1
       if P(i)>P(i-1)&&P(i)>P(i+1)
           estimate_theta(flag)=i;
           flag=flag+1;
       end
   else
       continue;
   end
end

%% 距离估计
P= zeros(target,length(R)); %波束方向图
estimate_distance=zeros(1,target); %RD-MUSIC算法距离维估计结果
for i=1:target
    for m = 1 : length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(estimate_theta(i)),R(m)); %导向矢量
         J=a'*En*En'*a; 
        P(i,m) =1/J;
    end
    figure(i+1);
    D=abs(P(i,:))/max(abs(P(i,:)));
    plot(R,10*log10(D));
    xlabel('距离/m');  ylabel('归一化幅度/dB'); 
    title('距离维估计结果');
    axis([0 Ru -40 0]);
    estimate_distance(i)=(find(D==max(D))*Ru-1)/1000;

end




