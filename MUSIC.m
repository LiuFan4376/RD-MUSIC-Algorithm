% 2D-MUSIC�㷨

clc; 
clear ; 
close ;

%% ---------------------------��������
j = sqrt(-1);
c = 3e8;
M= 8;             % ������Ԫ��
N = 6;             % ������Ԫ��
f0 = 2e10;          % �ο�Ƶ��
lamda0 = c/f0;
dt = lamda0/2;      % ������Ԫ���
dr = lamda0/2;      % ������Ԫ���
Dt=(0:M-1)*dt;      % ����������Ԫ�������
Dr=(0:N-1)*dr;      % ����������Ԫ������� 
delta_f= 1000;   % Ƶ�ʲ�����
Delta_f=(0:M-1)*delta_f; %��������Ƶƫ����
K = 300;   %������ 
target=2;% Ŀ����Ŀ
Ru=c/(2*delta_f);
theta=(-90:1:90)*pi/180; %�����Ƕ�����
R=linspace(0,Ru,Ru); %������������

%% ---------------------���������ź�
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %���������ź� 



%% ------------------------------------------- ������һ��Ŀ���ź�
SNR = 10;             
R1 = 1e5;          %delta_rΪ0
theta1 = 60/180*pi;                          %�Ƕ�
a1=steer_vector(f0,Delta_f,Dt,Dr,theta1,R1);  %Ŀ��һ�ĵ���ʸ��
sig = a1 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)��ɢ��ϵ��
sig1 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % ���������




%% ------------------------------------------- �����ڶ���Ŀ���ź�
SNR = 10;             
R2 = 0.75e5;          %delta_rΪ0
theta2 = 30/180*pi;                          %�Ƕ�
a2=steer_vector(f0,Delta_f,Dt,Dr,theta2,R2);  %Ŀ��һ�ĵ���ʸ��
sig = a2 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)��ɢ��ϵ��
sig2 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % ���������




%% --------------------------------------------�����źŴ���
data = sig1 +sig2+ noise;        %��������
% re = res_r*r_num;          % ������ƾ���
% g = kron(ones(N,1),exp(j*2*pi*2*delta_f/c*re*(0:M-1)'));        %����ʸ��
% data_comp = zeros(size(data));         
% for k = 1 : K
%     data_comp(:,k) = data(:,k).*g;             % ������Ľ�������
% end
% data_FDA_MIMO =data_comp ;   

Rxx=data*data'/K;

%% ����ֵ�ֽ�
[EV,D]=eig(Rxx);                   %����ֵ�ֽ�
EVA=diag(D)';                      %������ֵ����Խ�����ȡ��תΪһ��
[EVA,I]=sort(EVA);                 %������ֵ���� ��С����
EV=fliplr(EV(:,I));                % ��Ӧ����ʸ������
En=EV(:,target+1:M*N);                  %�����ӿռ�

%% �׷�����
P= zeros(length(theta),length(R)); %��������ͼ

 for n = 1 : length(theta)
    for m = 1 : length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(n),R(m)); %����ʸ��
         J=a'*En*En'*a;
        P(n,m) =1/J;
    end
 end
 
P=P';
% [x,y]=find(P==max(max(P)));
% esti_angle=y-91; %�㷨���ƽǶ�
% esti_distance=Ru/1000*(x-1); %�㷨���ƾ���
% disp(['target 1: angle =',esti_angle,' distance =',esti_distance]);


figure(1); 
imagesc(theta*180/pi,R,10*log(abs(P)/max(max(abs(P))))); 
xlabel('\theta^o'); ylabel('R/m'); 
axis tight; axis xy;
title('2D MUSIC�㷨���ƽ��');



