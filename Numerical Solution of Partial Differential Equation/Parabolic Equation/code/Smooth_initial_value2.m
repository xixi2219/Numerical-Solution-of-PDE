%% 显示格式解的图
mu1 = 1/4;
T0 = 1;
N = 256;
x = linspace(0,1,N+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_Explicit_Solver_1(N,mu1,U_0);% 数值解
u2 = exp(-pi^2*T0)*sin(pi*x(2:end-1));% 真解

figure(1)
subplot(1,2,1)
plot(x(2:end-1),u1,'-x','MarkerIndices',1:8:length(u1));
title(join(["N=",num2str(N,'%u'), "数值解"]));
subplot(1,2,2)
plot(x(2:end-1),u2,'-*','MarkerIndices',1:8:length(u2));
title(join(["N=",num2str(N,'%u'), "精确解"]));

%% 隐示格式解的图
mu1 = 1/4;
T0 = 1;
N = 256;
x = linspace(0,1,N+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_Implicit_Solver_1(N,mu1,U_0);% 数值解
u2 = exp(-pi^2*T0)*sin(pi*x(2:end-1));% 真解

figure(2)
subplot(1,2,1)
plot(x(2:end-1),u1,'-x','MarkerIndices',1:8:length(u1));
title(join(["N=",num2str(N,'%u'), "数值解"]));
subplot(1,2,2)
plot(x(2:end-1),u2,'-*','MarkerIndices',1:8:length(u2));
title(join(["N=",num2str(N,'%u'), "精确解"]));

%% CN格式解的图
mu1 = 1/4;
T0 = 1;
N = 256;
x = linspace(0,1,N+1);
U_0 = sin(pi*x(2:end-1));
u1 =  HeatEquation_CN_Solver_1(N,mu1,U_0);% 数值解
u2 = exp(-pi^2*T0)*sin(pi*x(2:end-1));% 真解

figure(3)
subplot(1,2,1)
plot(x(2:end-1),u1,'-x','MarkerIndices',1:8:length(u1));
title(join(["N=",num2str(N,'%u'), "数值解"]));
subplot(1,2,2)
plot(x(2:end-1),u2,'-*','MarkerIndices',1:8:length(u2));
title(join(["N=",num2str(N,'%u'), "精确解"]));