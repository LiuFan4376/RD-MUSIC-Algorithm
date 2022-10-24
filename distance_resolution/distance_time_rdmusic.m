clc;
close;
clear;
% % 时间已过 1.712997 秒。
% % 时间已过 0.727753 秒。
% % 时间已过 0.490123 秒。
% % 时间已过 0.378705 秒。
% % 时间已过 0.309404 秒。
% % 时间已过 0.263134 秒。
% % 时间已过 0.228112 秒。
% % 时间已过 0.201267 秒。
% % 时间已过 0.183322 秒。
% % 时间已过 0.167476 秒。
%% --------------------------参数设置
j = sqrt(-1);
c = 3e8;
M= 8;             %   发射阵元数
N = 8;              % 接收阵元数
f0 = 3e10;          % 参考频率
lamda0 = c/f0;
dt = lamda0/2;      % 发射阵元间距
dr = lamda0/2;      % 接收阵元间距
Dt=(0:M-1)*dt;      % 发射阵列阵元间距设置
Dr=(0:N-1)*dr;      % 接收阵列阵元间距设置 
delta_f= 1000;   % 频率步进量
Delta_f=(0:M-1)*delta_f; %发射阵列频偏设置
K = 500;   %快拍数
target_num=1;% 目标数目
Ru=c/(2*delta_f);

target_angle=[0];  %目标角度
target_distance=[5e3];   % 目标距离

distance_resolution=1:1:10; %距离维分辨率
theta=(-90:1:90)*pi/180;
 
% target_angle=[0.9]; 
% target_distance=[0.5001e5];   

snr=10;   

%% ---------------------------------------接收信号处理

%%  噪声设置
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %浜х????澹颁俊?? 


%% 目标信号设置
data=noise;  
for i=1:target_num          
    a=steer_vector(f0,Delta_f,Dt,Dr,target_angle(i)*pi/180,target_distance(i));  %导向矢量
    sig = a * sqrt(10^(snr/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)
    sig= sqrt(10^(snr/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 
    data=data+sig;
end

for l=distance_resolution

    tic;
 R=1:l:2e4;
 
%% 矩阵特征值分解
Rxx=data*data'/K;                
 G=Rxx(:,1:target_num);
H=Rxx(:,target_num+1:M*N);
P=inv(G'*G)*G'*H;
I=-diag(ones(1,M*N-target_num));
En=[P;I];
En=orth(En);
En = En./repmat(sqrt(sum(En.^2,1)),size(En,1),1);


%% ---------------------------------------------RD-MUSIC

%% 角度估计
P= zeros(1,length(theta)); %角度维波束图
 for n = 1:length(theta)
         
        d = exp(j*2*pi*f0/c*Dt'*sin(theta(n)));     %  发射角度导向矢量
         b = exp(j*2*pi*f0/c*Dr'*sin(theta(n)));     %  接收导向矢量
          W=kron(b,diag(d))'*(En*En')*kron(b,diag(d));
          W1=W(1,1);
          W2=W(1,2:M);
          W4=W(2:M,2:M);
          J=W1-W2/W4*W2';
        P(n) =1/J;
   
 end
 P=abs(P)/max(abs(P));
estimate_theta=zeros(1,target_num); %RD-MUSIC目标角度估计结果
% flag=1;
% for i=1:length(P)
%     if flag==target_num+1
%         break;
%     end
%     if P(i)<=0.1
%         continue;
%     end
%  
%     if i==1&&P(i)>P(i+1)
%         estimate_theta(flag)=i;
%         flag=flag+1;
%     end
%     if i==length(P)&&P(i)>P(i-1)
%         estimate_theta(flag)=i;
%         flag=flag+1;
%     end
%    if i<length(P)&&i>1
%        if P(i)>P(i-1)&&P(i)>P(i+1)
%            estimate_theta(flag)=i;
%            flag=flag+1;
%        end
%    else
%        continue;
%    end
% end

    estimate_theta(1)=find(P==max(P));
% 距离估计
 P= zeros(target_num,length(R)); %娉距离维波束方向图
estimate_distance=zeros(1,target_num); %距离维估计结果
for i=1:target_num
    for m =1:length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(estimate_theta(i)),R(m)); %联合导向矢量
         J=a'*(En*En')*a; 
        P(i,m) =1/J;
    end
    D=abs(P(i,:))/max(abs(P(i,:)));
    estimate_distance(i)=find(D==max(D));

end
  toc;
end
