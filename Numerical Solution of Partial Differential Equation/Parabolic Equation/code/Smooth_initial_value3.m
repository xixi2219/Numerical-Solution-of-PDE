
%% 显式格式不稳定情况
mu1 = 0.502;
T0 = 1;
N = [32 64 128 256];
for i=1:4
h = 1/N(i);
Delta_t = h^2*mu1;
NumT = T0/ Delta_t;
NumT = fix(NumT);
T1 = NumT*Delta_t;
x = linspace(0,1,N(i)+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_Explicit_Solver_1(N(i),mu1,U_0);% 数值解
u2 = exp(-pi^2*T1)*sin(pi*x(2:end-1));% 真解
figure(1)
subplot(2,2,i)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('数值结果','真解');
title(join(["N=",num2str(N(i),'%u')]));
end

%% 隐式格式mu很大情况
mu1 = 10;
T0 = 1;
N = 256;
h = 1/N;
Delta_t = h^2*mu1;
NumT = T0/ Delta_t;
NumT = fix(NumT);
T1 = NumT*Delta_t;
x = linspace(0,1,N+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_Implicit_Solver_1(N,mu1,U_0);% 数值解
u2 = exp(-pi^2*T1)*sin(pi*x(2:end-1));% 真解
figure(2)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('数值结果','真解');
title(join(["N=",num2str(N,'%u')]));

%% CN格式mu很大情况
mu1 = 100;
T0 = 1;
N = [32 64 128 256];
for i=1:4
h = 1/N(i);
Delta_t = h^2*mu1;
NumT = T0/ Delta_t;
NumT = fix(NumT);
T1 = NumT*Delta_t;
x = linspace(0,1,N(i)+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_CN_Solver_1(N(i),mu1,U_0);% 数值解
u2 = exp(-pi^2*T1)*sin(pi*x(2:end-1));% 真解
figure(3)
subplot(2,2,i)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('数值结果','真解');
title(join(["N=",num2str(N(i),'%u')]));
end
