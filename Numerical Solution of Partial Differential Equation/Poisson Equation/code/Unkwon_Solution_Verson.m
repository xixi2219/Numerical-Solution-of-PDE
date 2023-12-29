%% 无精确解的情况
e=1:8;
N=[4,8,16,32,64,128,256,512];
t=1;
for i=1:8
    h = 1/N(i);
    x = 0:h:1; y = 0:h:1;
    F =  Exponent_Func(x , y);
    M =  Dirichlet_partial(x , y);
   if i==1
        [U1h] = PossionEquation1(N(i),F,M);
    else
        [U2h] = PossionEquation1(N(i),F,M);
        middleU=U2h(2:2:N(i)-1,2:2:N(i)-1);
        e2(i) = max(max(abs(U1h-middleU)));
        U1h = U2h;
   end
   if ismember(N(i),[32,128,512])
        [X,Y] = meshgrid(x,y);%格点坐标
        figure(1)
        subplot(1,3,t)
        mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),U1h)
        view(3)
        title(join(["N=",num2str(N(i),'%u'),"数值解"]))
        t=t+1;
   end
end

%% 误差界
figure

plot(log(1./N(2:end)),log(e2(2:end)),'-s')
[a,b]=polyfit(log(1./N(2:end)),log(e2(2:end)),1);

legend('$\ln||U_{h/2}-U_{h}|| = 1.9674\ln h+C$','Interpreter','latex');

xlabel('$\ln h$','Interpreter','latex')
grid on
ylabel('$\ln||U_{h/2}-U_{h}||$','Interpreter','latex')





%% 子函数
function [f] =  Exponent_Func(x , y)
% -f在网格点上的值
    [X,Y] = meshgrid(x,y);
    f = -(X.*X+Y.*Y).*exp(X+Y);
end

function [M] =  Dirichlet_partial(x , y)
% 边界值
    [X,Y] = meshgrid(x,y);
    M = X.*X+Y.*Y;
    M(2:end-1,2:end-1) = 0;
end

