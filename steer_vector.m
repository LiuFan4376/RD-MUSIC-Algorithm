function a = steer_vector(f0,Delta_f,Dt,Dr,theta,R)
%steer_vector ����FDA-MIMO�״�ĵ���ʸ��
%   f0Ϊ�ز�����Ƶ��,DtΪ�������м�����ã�DrΪ���������м�����ã�Delta_fΪ��������������ԪƵƫ���ã�theta��RΪλ����Ϣ


%% -------------------------------------�״��������
j = sqrt(-1);
c = 3e8;

%% -------------------------------------����ʸ��
a_t_r= exp(-j*2*pi*2*Delta_f'*R/c);            %  �������о��뵼��ʸ��
a_t_theta = exp(j*2*pi*f0/c*Dt'*sin(theta));        %  �������нǶȵ���ʸ��
b = exp(j*2*pi*f0/c*Dr'*sin(theta));        %  �������е���ʸ��
a=kron(b,a_t_r.*a_t_theta);                 %FDA-MIMOϵͳ����ʸ��



end

