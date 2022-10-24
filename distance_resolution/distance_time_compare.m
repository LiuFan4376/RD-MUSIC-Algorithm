clc;
close;
clear;

distance_resolution=1:1:10; %距离维分辨率
distance_time_2dmusic=[131.723541  33.498801  14.670724  8.249863  5.492221  4.098453  2.711773  2.066746  1.693376 1.367774  ];
distance_time_rdmusic=[1.712997 0.727753  0.490123 0.378705  0.309404  0.263134  0.228112  0.201267  0.183322 0.167476 ];
%% 画图   
 figure(1);
plot(distance_resolution,distance_time_2dmusic,'s-',distance_resolution,distance_time_rdmusic,'o-','LineWidth',1.5);
xlabel('距离维划分精度/m'); ylabel('时间/s'); 
% axis([-90,90,0,1]);
legend('2D-MUISC','RD-MUSIC');