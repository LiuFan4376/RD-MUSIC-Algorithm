function a = steer_vector(f0,Delta_f,Dt,Dr,theta,R)
%steer_vector 产生FDA-MIMO雷达的导向矢量
%   f0为载波中心频率,Dt为发射阵列间距设置，Dr为接收阵阵列间距设置，Delta_f为发射阵列相邻阵元频偏设置，theta和R为位置信息


%% -------------------------------------雷达参数设置
j = sqrt(-1);
c = 3e8;

%% -------------------------------------导向矢量
a_t_r= exp(-j*2*pi*2*Delta_f'*R/c);            %  发射阵列距离导向矢量
a_t_theta = exp(j*2*pi*f0/c*Dt'*sin(theta));        %  发射阵列角度导向矢量
b = exp(j*2*pi*f0/c*Dr'*sin(theta));        %  接收阵列导向矢量
a=kron(b,a_t_r.*a_t_theta);                 %FDA-MIMO系统导向矢量



end

