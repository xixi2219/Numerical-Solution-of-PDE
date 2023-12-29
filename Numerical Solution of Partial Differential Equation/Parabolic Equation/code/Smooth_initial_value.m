%% Space Error Analysis
mu = 1/4;
N = [32,64,128,256,512,1024];

%Explicit
for i =1:length(N)
   [Error_Explicit1(i),Error_Explicit2(i)] = HeatEquation_Explicit_Solver(N(i),mu);
end
%收敛阶
Order_Explicit1 = (log(Error_Explicit1(2:end))-log(Error_Explicit1(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));
Order_Explicit2 = (log(Error_Explicit2(2:end))-log(Error_Explicit2(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));

%Implicit
  for i=1:length(N)
      [Error_Implicit1(i),Error_Implicit2(i)] = HeatEquation_Implicit_Solver(N(i),mu);
  end
  Order_Implicit1 = (log(Error_Implicit1(2:end))-log(Error_Implicit1(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));
  Order_Implicit2 = (log(Error_Implicit2(2:end))-log(Error_Implicit2(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));
  
% CN
  for i=1:length(N)
      [Error_CN1(i),Error_CN2(i)] =HeatEquation_CN_Solver(N(i),mu);
  end
  Order_CN1 = (log(Error_CN1(2:end))-log(Error_CN1(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));
  Order_CN2 = (log(Error_CN2(2:end))-log(Error_CN2(1:end-1)))./(log(N(1:end-1))-log(N(2:end)));

%显格式图
figure(1)
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
plot(log(1./N),log(Error_Explicit1),'black--o')
[a_ex1,b_ex1]=polyfit(log(1./N),log(Error_Explicit1),1);
legend('$\ln ||e||_{\infty}=1.9999\ln h + C$','Interpreter','latex')



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

plot(log(1./N),log(Error_Explicit2),'black--o')
[a_ex2,b_ex2]=polyfit(log(1./N),log(Error_Explicit2),1);
legend('$\ln ||e||_{2}=1.9999\ln h + C$','Interpreter','latex')

%隐格式图
figure(2)
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
plot(log(1./N),log(Error_Implicit1),'black--o')
[a_im1,b_im1]=polyfit(log(1./N),log(Error_Implicit1),1);
legend('$\ln ||e||_{\infty}=2.0021\ln h + C$','Interpreter','latex')


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
plot(log(1./N),log(Error_Implicit2),'black--o')
[a_im2,b_im2]=polyfit(log(1./N),log(Error_Implicit2),1);
legend('$\ln ||e||_{2}=2.0021\ln h + C$','Interpreter','latex')

%CN格式图
figure(3)
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
plot(log(1./N),log(Error_CN1),'black--o')
[a_cn1,b_cn1]=polyfit(log(1./N),log(Error_CN1),1);
legend('$\ln ||e||_{\infty}=2.0009\ln h + C$','Interpreter','latex')


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
plot(log(1./N),log(Error_CN2),'black--o')
[a_cn2,b_cn2]=polyfit(log(1./N),log(Error_CN2),1);
legend('$\ln ||e||_{2}=2.0009\ln h + C$','Interpreter','latex')