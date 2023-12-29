% error analysis of the three methods
T = 1;
meshrate = 0.8;
N = [64,128,256,512,1024];
Error_infty = zeros(1,length(N));
Error_L2 = zeros(1,length(N));
for i = 1:length(N)
x = linspace(0,1,N(i)+1);
h = 1/N(i);
u0 = sin(2*pi.*x);% Initial Value
[u_T,T_max] = HyperbolicEquationSolver(u0,h,meshrate,T,3);%三种格式都可计算
u_true = sin(2*pi.*(x-T_max));% 真解
Error_infty(i) = max(abs(u_true-u_T));
Error_L2(i) = sqrt( 1/(N(i)/2+1/2) * sum(abs(u_true-u_T).^2));
end
% 计算收敛阶
Order_infty = (log(Error_infty(2:end))-log(Error_infty(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));
Order_L2 = (log(Error_L2(2:end))-log(Error_L2(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));

% 收敛阶图像
figure(1)
sgtitle('$\nu=0.8$','Interpreter','latex')

subplot(1,2,1)
hold on
grid on

ax = gca;
ax.GridLineStyle=":";
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = ":";
ax.Box = 'on';
ax.LineWidth = 0.75;
xlabel('$\ln h$','Interpreter','latex')
ylabel('$\ln||e||_{\infty}$','Interpreter','latex')
plot(log(1./N),log(Error_infty),'black--o')
[a_ex1,b_ex1]=polyfit(log(1./N),log(Error_infty),1);
legend('$\ln ||e||_{\infty}=1.9998\ln h + C$','Interpreter','latex')

subplot(1,2,2)
hold on
grid on

ax = gca;
ax.GridLineStyle=":";
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = ":";
ax.Box = 'on';
ax.LineWidth = 0.75;
xlabel('$\ln h$','Interpreter','latex')
ylabel('$\ln||e||_{2}$','Interpreter','latex')

plot(log(1./N),log(Error_L2),'black--o')
[a_ex2,b_ex2]=polyfit(log(1./N),log(Error_L2),1);
legend('$\ln ||e||_{2}=2.0023\ln h + C$','Interpreter','latex')