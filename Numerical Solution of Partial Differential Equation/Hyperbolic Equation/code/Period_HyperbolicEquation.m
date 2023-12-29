%% 稳定情况
N = 512;
h = 1/N;
meshrate = 1.5;
x = linspace(0,1,N+1);
u0 = sin(2*pi.*x);% Initial Value

figure(1)
subplot(1,3,1)
plot(x,u0);
title('$t=0$','Interpreter','latex');
legend('真解');

subplot(1,3,2)
[u_T1,T_max1] = HyperbolicEquationSolver(u0,h,meshrate,0.5,3);%三种格式都可计算
u_true1 = sin(2*pi.*(x-T_max1));% 真解
plot(x,u_true1,x,u_T1,'-.x','MarkerIndices',1:8:length(u_T1));
title('$t=0.5$','Interpreter','latex');
legend('真解','数值解');

subplot(1,3,3)
[u_T2,T_max2] = HyperbolicEquationSolver(u0,h,meshrate,1,3);%三种格式都可计算
u_true2 = sin(2*pi.*(x-T_max2));% 真解
plot(x,u_true2,x,u_T2,'-.x','MarkerIndices',1:8:length(u_T2));
title('$t=1$','Interpreter','latex');
legend('真解','数值解');

%% 不稳定情况
N = 512;
h = 1/N;
meshrate = 2.5;
x = linspace(0,1,N+1);
u0 = sin(2*pi.*x);% Initial Value

figure(2)
subplot(1,3,1)
plot(x,u0);
title('$t=0$','Interpreter','latex');
legend('真解');

subplot(1,3,2)
[u_T3,T_max3] = HyperbolicEquationSolver(u0,h,meshrate,0.5,3);%三种格式都可计算
u_true3 = sin(2*pi.*(x-T_max3));% 真解
plot(x,u_true3,x,u_T3,'-.x','MarkerIndices',1:8:length(u_T3));
title('$t=0.5$','Interpreter','latex');
legend('真解','数值解');

subplot(1,3,3)
[u_T4,T_max4] = HyperbolicEquationSolver(u0,h,meshrate,1,3);%三种格式都可计算
u_true4 = sin(2*pi.*(x-T_max4));% 真解
plot(x,u_true4,x,u_T4,'-.x','MarkerIndices',1:8:length(u_T4));
title('$t=1$','Interpreter','latex');
legend('真解','数值解');

%% 色散与耗散
N = 128;
h = 1/N;
meshrate = 0.8;
x = linspace(0,1,N+1);
u0 = sin(2*pi.*x);% Initial Value

figure(3)
subplot(1,3,1)
[u_T5,T_max5] = HyperbolicEquationSolver(u0,h,meshrate,100,3);%三种格式都可计算
u_true5 = sin(2*pi.*(x-T_max5));% 真解
plot(x,u_true5,x,u_T5,'-.x','MarkerIndices',1:8:length(u_T5));
title('$t=100$','Interpreter','latex');
legend('真解','数值解');

subplot(1,3,2)
[u_T6,T_max6] = HyperbolicEquationSolver(u0,h,meshrate,150,3);%三种格式都可计算
u_true6 = sin(2*pi.*(x-T_max6));% 真解
plot(x,u_true6,x,u_T6,'-.x','MarkerIndices',1:8:length(u_T6));
title('$t=150$','Interpreter','latex');
legend('真解','数值解');

subplot(1,3,3)
[u_T7,T_max7] = HyperbolicEquationSolver(u0,h,meshrate,250,3);%三种格式都可计算
u_true7 = sin(2*pi.*(x-T_max7));% 真解
plot(x,u_true7,x,u_T7,'-.x','MarkerIndices',1:8:length(u_T7));
title('$t=250$','Interpreter','latex');
legend('真解','数值解');