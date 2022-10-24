%% RD-MUSIC�㷨
clc;
close;
clear;

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
target=3;% Ŀ����Ŀ
Ru=c/(2*delta_f);
theta=(-90:0.1:90)*pi/180; %�����Ƕ�����
R=linspace(0,Ru,1000); %������������

%% ---------------------���������ź�
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %���������ź� 

%% ------------------------------------------- ������һ��Ŀ���ź�
SNR = 10;             
R1 = 10e4;          %delta_rΪ0
theta1 = -20/180*pi;                          %�Ƕ�
a1=steer_vector(f0,Delta_f,Dt,Dr,theta1,R1);  %Ŀ��һ�ĵ���ʸ��
sig = a1 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)��ɢ��ϵ��
sig1 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % ���������

%% ------------------------------------------- �����ڶ���Ŀ���ź�
SNR = 10;             
R2 = 7.5e4;          %delta_rΪ0
theta2 = -10/180*pi;                          %�Ƕ�
a2=steer_vector(f0,Delta_f,Dt,Dr,theta2,R2);  %Ŀ��һ�ĵ���ʸ��
sig = a2 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)��ɢ��ϵ��
sig2 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % ���������
SNR = 10;
%% ------------------------------------------- ����������Ŀ���ź�
SNR = 10;
R3 = 8.5e4;          %delta_rΪ0
theta3 = 20/180*pi;                          %�Ƕ�
a3=steer_vector(f0,Delta_f,Dt,Dr,theta3,R3);  %Ŀ��һ�ĵ���ʸ��
sig = a3 * sqrt(10^(SNR/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)��ɢ��ϵ��
sig3 = sqrt(10^(SNR/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % ���������

%% --------------------------------------------�����źŴ���
data = sig1 +sig2+sig3+noise;        %��������  
 Rxx=data*data'/K;                %Э�������

%% ����ֵ�ֽ�
% [EV,D]=eig(Rxx);                   %����ֵ�ֽ�
% EVA=diag(D)';                      %������ֵ����Խ�����ȡ��תΪһ��
% [EVA,I]=sort(EVA);                 %������ֵ���� ��С����
% EV=fliplr(EV(:,I));                % ��Ӧ����ʸ������
% En=EV(:,target+1:M*N);                  %�����ӿռ�
G=Rxx(:,1:target);
H=Rxx(:,target+1:M*N);
P=inv(G'*G)*G'*H;
I=-diag(ones(1,M*N-target));
En=[P;I];
En=orth(En);
En = En./repmat(sqrt(sum(En.^2,1)),size(En,1),1);


%% -----------------------------------------------------�Ƕȹ���
P= zeros(1,length(theta)); %��������ͼ
 for n = 1 : length(theta)
         
        d = exp(j*2*pi*f0/c*Dt'*sin(theta(n)));        %  �������нǶȵ���ʸ��
         b = exp(j*2*pi*f0/c*Dr'*sin(theta(n)));        %  �������е���ʸ�� 
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
xlabel('�Ƕ�/o');  ylabel('��һ������/dB'); 

axis([-90.1 90.1 -35 0]);
%% �����Ƕȹ���ֵ
estimate_theta=zeros(1,target); %RD-MUSIC�㷨�Ƕ�ά���ƽ��
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

%% �������
P= zeros(target,length(R)); %��������ͼ
estimate_distance=zeros(1,target); %RD-MUSIC�㷨����ά���ƽ��
for i=1:target
    for m = 1 : length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(estimate_theta(i)),R(m)); %����ʸ��
         J=a'*En*En'*a; 
        P(i,m) =1/J;
    end
    figure(i+1);
    D=abs(P(i,:))/max(abs(P(i,:)));
    plot(R,10*log10(D));
    xlabel('����/m');  ylabel('��һ������/dB'); 
    title('����ά���ƽ��');
    axis([0 Ru -40 0]);
    estimate_distance(i)=(find(D==max(D))*Ru-1)/1000;

end




