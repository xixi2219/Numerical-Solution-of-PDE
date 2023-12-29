%% 方波初值与三角波初值情况
N = 64*310;
h = 310/N;
meshrate = 0.8;
x = linspace(0,310,N+1);
% u0 = squarewave(x);% Initial Value
u0 = triwave(x);

figure(1)
subplot(2,3,1)
plot(x(1:64*5+1),u0(1:64*5+1));
title('$t=0$','Interpreter','latex');
legend('真解');

subplot(2,3,2)
[u_T1,T_max1] = HyperbolicEquationSolver(u0,h,meshrate,3,3);%三种格式都可计算
% u_true1 = squarewave(x-T_max1);% 真解
u_true1 = triwave(x-T_max1);
plot(x(64*1+1:64*6+1),u_true1(64*1+1:64*6+1),x(64*1+1:64*6+1),u_T1(64*1+1:64*6+1),'-.x','MarkerIndices',1:8:length(u_T1(64*1+1:64*6+1)));
title('$t=3$','Interpreter','latex');
legend('真解','数值解');

subplot(2,3,3)
[u_T2,T_max2] = HyperbolicEquationSolver(u0,h,meshrate,5,3);%三种格式都可计算
% u_true2 = squarewave(x-T_max2);% 真解
u_true2 = triwave(x-T_max2);
plot(x(64*3+1:64*8+1),u_true2(64*3+1:64*8+1),x(64*3+1:64*8+1),u_T2(64*3+1:64*8+1),'-.x','MarkerIndices',1:8:length(u_T2(64*3+1:64*8+1)));
title('$t=5$','Interpreter','latex');
legend('真解','数值解');

subplot(2,3,4)
[u_T3,T_max3] = HyperbolicEquationSolver(u0,h,meshrate,100,3);%三种格式都可计算
% u_true3 = squarewave(x-T_max3);% 真解
u_true3 = triwave(x-T_max3);
plot(x(64*98+1:64*103+1),u_true3((64*98+1:64*103+1)),x(64*98+1:64*103+1),u_T3(64*98+1:64*103+1),'-.x','MarkerIndices',1:8:length(u_T3(64*98+1:64*103+1)));
title('$t=100$','Interpreter','latex');
legend('真解','数值解');

subplot(2,3,5)
[u_T4,T_max4] = HyperbolicEquationSolver(u0,h,meshrate,200,3);%三种格式都可计算
% u_true4 = squarewave(x-T_max4);% 真解
u_true4 = triwave(x-T_max4);
plot(x(64*198+1:64*203+1),u_true4(64*198+1:64*203+1),x(64*198+1:64*203+1),u_T4(64*198+1:64*203+1),'-.x','MarkerIndices',1:8:length(u_T4(64*198+1:64*203+1)));
title('$t=200$','Interpreter','latex');
legend('真解','数值解');

subplot(2,3,6)
[u_T5,T_max5] = HyperbolicEquationSolver(u0,h,meshrate,300,3);%三种格式都可计算
% u_true5 = squarewave(x-T_max5);% 真解
u_true5 = triwave(x-T_max5);
plot(x(64*298+1:64*303+1),u_true5(64*298+1:64*303+1),x(64*298+1:64*303+1),u_T5(64*298+1:64*303+1),'-.x','MarkerIndices',1:8:length(u_T5(64*298+1:64*303+1)));
title('$t=300$','Interpreter','latex');
legend('真解','数值解');






%% 初值函数
% 方波初值
function u0 = squarewave(x)
u0 = zeros(1,length(x));
for i = 1:length(x)
    if x(i) == 0
        u0(i) = 0;
    elseif (x(i)>0)&&(x(i)<=1)
        u0(i) = 1;
    else
        u0(i) = 0;
    end
end
end
% 三角波初值
function u0 = triwave(x)
u0 = zeros(1,length(x));
for i = 1:length(x)
    if x(i) == 0
        u0(i) = 0;
    elseif (x(i)>0)&&(x(i)<=1)
        u0(i) = x(i);
    elseif (x(i)>1)&&(x(i)<=2)
        u0(i) = 2-x(i);
    else 
        u0(i) = 0;
    end
end
end