
%% CN格式非光滑初值不同mu的对照分析
mu1 = 1/4;
mu2 = 100;
N = [64 128 256 512];
for i = 1:length(N)
x = linspace(0,1,N(i)+1);
U_0 = initial_condition(x(2:end-1));
u1 =  HeatEquation_CN_Solver_1(N(i),mu1,U_0);
u2 = HeatEquation_CN_Solver_1(N(i),mu2,U_0);

figure(1)
subplot(2,2,i)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('真解','数值结果');
title(join(["N=",num2str(N(i),'%u')]));
end

%% 显格式非光滑初值不同mu的对照分析
mu1 = 1/4;
mu2 = 0.502;
N = [32 64 128 256];
for i = 1:length(N)
x = linspace(0,1,N(i)+1);
U_0 = initial_condition(x(2:end-1));
u1 =  HeatEquation_Explicit_Solver_1(N(i),mu1,U_0);
u2 = HeatEquation_Explicit_Solver_1(N(i),mu2,U_0);

figure(2)
subplot(2,2,i)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('真解','数值结果');
title(join(["N=",num2str(N(i),'%u')]));
end

%% 初值
figure(3)
plot(x(2:end-1),U_0);
title('非光滑初值');

%% 隐格式非光滑初值不同mu的对照分析
mu1 = 1/4;
mu2 = 10;
N = [32 64 128 256];
for i = 1:length(N)
x = linspace(0,1,N(i)+1);
U_0 = initial_condition(x(2:end-1));
u1 =  HeatEquation_Implicit_Solver_1(N(i),mu1,U_0);
u2 = HeatEquation_Implicit_Solver_1(N(i),mu2,U_0);

figure(4)
subplot(2,2,i)
hold on
grid on
plot(x(2:end-1),u1,'-');
plot(x(2:end-1),u2,'-.');
legend('真解','数值结果');
title(join(["N=",num2str(N(i),'%u')]));
end







%% 非光滑初值
function U_0 = initial_condition(x)
for i = 1:length(x)
    if (x(i)<=(1/4))&&(x(i)>=0)
        U_0(i) = 4.*x(i);
    elseif (x(i)>1/4)&&(x(i)<=1/2)
        U_0(i) = -4.*x(i)+2;
    elseif (x(i)>1/2)&&(x(i)<=3/4)
        U_0(i) = 4.*x(i)-2;
    elseif (x(i)>3/4)&&(x(i)<=1)
        U_0(i) = -4.*x(i)+4;  
    else
        error("x的值不符合(0,1)的要求");
    end
end
end