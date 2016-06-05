clc,close all,


load ss1.mat
load ss2.mat
load ss3.mat
load ss4.mat

plot(ss1(1,:),ss1(2,:),'^-',ss2(1,:),ss2(2,:),'*-',ss3(1,:),ss3(2,:),'o-',...
    ss4(1,:),ss4(2,:),'x-')

legend('Re=1','Re=100','Re=1000','Re=2000')

grid on
xlabel('Node number from left to right for the centerline')
ylabel('v_y')