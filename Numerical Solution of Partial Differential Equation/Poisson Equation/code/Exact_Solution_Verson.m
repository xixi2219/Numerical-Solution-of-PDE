%% 有精确解的情况
e=1:8;
N=[4,8,16,32,64,128,256,512];
t=1;
for i=1:8
    h = 1/N(i);
    x = 0:h:1; y = 0:h:1;
    F = Trigonometric_Func(x , y);
    M =  Dirichlet_partial(x , y);
    u = Exact_Solution(x, y);
   [U,e(i)]= PossionEquation1(N(i),F,M,u);
   if ismember(N(i),[32,128,512])
        [X,Y] = meshgrid(x,y);%格点坐标
        figure(1)
        subplot(1,3,t)
        mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),u(2:end-1,2:end-1))
        view(3)
        title(join(["N=",num2str(N(i),'%u'), "精确解"]))
        
        figure(2)
        subplot(1,3,t)
        mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),U)
        view(3)
        title(join(["N=",num2str(N(i),'%u'),"数值解"]))
        t=t+1;
   end
end
%估计收敛阶log_{h2/h1}^{e2/e1}
order=1:7;
for i=1:7
    order(i) = log(e(i+1)/e(i))/log(1/2);
end

% 误差界
figure

plot(log(1./N),log(e),'-s')
[a,b]=polyfit(log(1./N),log(e),1);
legend('$\ln e =1.995\ln h + C$','Interpreter','latex')
xlabel('$\ln h$','Interpreter','latex')
grid on
ylabel('$\ln e$','Interpreter','latex')
    

%% 子函数
function [f] =  Trigonometric_Func(x , y)
% -f在网格点上的值
    [X,Y] = meshgrid(x,y);
    f = -2*(pi^2).*sin(pi*X).*cos(pi*Y);
end

function [M] =  Dirichlet_partial(x , y)
% 边界值
    [X,Y] = meshgrid(x,y);
    M = sin(pi*X).*cos(pi*Y);
    M(2:end-1,2:end-1) = 0;
end


function [u] = Exact_Solution(x, y)
% 精确解
    [X,Y] = meshgrid(x,y);
    u = sin(pi*X).*cos(pi*Y);
end

