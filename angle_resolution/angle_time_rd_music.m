clc;
close;
clear;

%% --------------------------��������
j = sqrt(-1);
c = 3e8;
M= 8;             %   ������Ԫ��
N = 8;              % ������Ԫ��
f0 = 3e10;          % �ο�Ƶ��
lamda0 = c/f0;
dt = lamda0/2;      % ������Ԫ���
dr = lamda0/2;      % ������Ԫ���
Dt=(0:M-1)*dt;      % ����������Ԫ�������
Dr=(0:N-1)*dr;      % ����������Ԫ������� 
delta_f= 1000;   % Ƶ�ʲ�����
Delta_f=(0:M-1)*delta_f; %��������Ƶƫ����
K = 500;   %������
target_num=1;% Ŀ����Ŀ
Ru=c/(2*delta_f);

target_angle=[0];  %Ŀ��Ƕ�
target_distance=[5e3];   % Ŀ�����

angle_resolution=0.1:0.1:1; %�Ƕ�ά�ֱ���

 
% target_angle=[0.9]; 
% target_distance=[0.5001e5];   

snr=10; 



   

%% ---------------------------------------�����źŴ���

%%  ��������
noise = 1/sqrt(2)*(randn(N*M,K)+j*randn(N*M,K));
noise = 1/sqrt(trace(noise*noise'/K)/(M*N)) * noise; %产�????声信?? 


%% Ŀ���ź�����
data=noise;  
for i=1:target_num          
    a=steer_vector(f0,Delta_f,Dt,Dr,target_angle(i)*pi/180,target_distance(i));  %����ʸ��
    sig = a * sqrt(10^(snr/10))*randn(1,K);      %sqrt(10^(SNR/10))*randn(1,K)
    sig= sqrt(10^(snr/10)/(trace(sig*sig'/K)/trace(noise*noise'/K))) * sig;    % 
    data=data+sig;
end

for l=angle_resolution

    tic;
 R=1:1e4;
theta=(-90:l:90)*pi/180; 
%% ��������ֵ�ֽ�
Rxx=data*data'/K;                
 G=Rxx(:,1:target_num);
H=Rxx(:,target_num+1:M*N);
P=inv(G'*G)*G'*H;
I=-diag(ones(1,M*N-target_num));
En=[P;I];
En=orth(En);
En = En./repmat(sqrt(sum(En.^2,1)),size(En,1),1);


%% ---------------------------------------------RD-MUSIC

%% �Ƕȹ���
P= zeros(1,length(theta)); %�Ƕ�ά����ͼ
 for n = 1:length(theta)
         
        d = exp(j*2*pi*f0/c*Dt'*sin(theta(n)));     %  ����Ƕȵ���ʸ��
         b = exp(j*2*pi*f0/c*Dr'*sin(theta(n)));     %  ���յ���ʸ��
          W=kron(b,diag(d))'*(En*En')*kron(b,diag(d));
          W1=W(1,1);
          W2=W(1,2:M);
          W4=W(2:M,2:M);
          J=W1-W2/W4*W2';
        P(n) =1/J;
   
 end
 P=abs(P)/max(abs(P));
estimate_theta=zeros(1,target_num); %RD-MUSICĿ��Ƕȹ��ƽ��


    estimate_theta(1)=find(P==max(P));
% �������
 P= zeros(target_num,length(R)); %泾���ά��������ͼ
estimate_distance=zeros(1,target_num); %����ά���ƽ��
for i=1:target_num
    for m =1:length(R)
         a= steer_vector(f0,Delta_f,Dt,Dr,theta(estimate_theta(i)),R(m)); %���ϵ���ʸ��
         J=a'*(En*En')*a; 
        P(i,m) =1/J;
    end
    D=abs(P(i,:))/max(abs(P(i,:)));
    estimate_distance(i)=find(D==max(D));

end
  toc;
end



