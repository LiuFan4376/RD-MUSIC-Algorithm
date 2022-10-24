clc; 
clear ; 
close ;
% % ʱ���ѹ� 131.723541 �롣
% % ʱ���ѹ� 33.498801 �롣
% % ʱ���ѹ� 14.670724 �롣
% % ʱ���ѹ� 8.249863 �롣
% % ʱ���ѹ� 5.492221 �롣
% % ʱ���ѹ� 4.098453 �롣
% % ʱ���ѹ� 2.711773 �롣
% % ʱ���ѹ� 2.066746 �롣
% % ʱ���ѹ� 1.693376 �롣
% % ʱ���ѹ� 1.367774 ��
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
K = 500;   %������ 
target_num=1;% Ŀ����Ŀ
target_angle=[0];  %Ŀ��Ƕ�
target_distance=[0.5e3];   % Ŀ�����

distance_resolution=1:1:10; %����ά�ֱ���
theta=(-90:1:90)*pi/180;


%%  ��������
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; 


%% Ŀ���ź�����
snr=10;
data=noise;  
for i=1:target_num          
    a=steer_vector(f0,Delta_f,Dt,Dr,target_angle(i)*pi/180,target_distance(i));  %����ʸ��
    sig = a * sqrt(10^(snr/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)
    sig= sqrt(10^(snr/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 
    data=data+sig;
end


%% --------------------------------------------�����źŴ���

for l=distance_resolution
 tic;
  R=1:l:2e4;
Rxx=data*data'/K;

%% ����ֵ�ֽ�
[EV,D]=eig(Rxx);                   %����ֵ�ֽ�
EVA=diag(D)';                      %������ֵ����Խ�����ȡ��תΪһ��
[EVA,I]=sort(EVA);                 %������ֵ���� ��С����
EV=fliplr(EV(:,I));                % ��Ӧ����ʸ������
En=EV(:,target_num+1:M*N);                  %�����ӿռ�



theta=(-90:l:90)*pi/180; 
    
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
[x,y]=find(P==max(max(P)));
% esti_angle=y-91; %�㷨���ƽǶ�
% esti_distance=Ru/1000*(x-1); %�㷨���ƾ���
% disp(['target 1: angle =',esti_angle,' distance =',esti_distance]);

 toc;
end
% 
% figure(1); 
% imagesc(theta*180/pi,R,10*log(abs(P)/max(max(abs(P))))); 
% xlabel('\theta^o'); ylabel('R/m'); 
% axis tight; axis xy;
% title('2D MUSIC�㷨���ƽ��');
