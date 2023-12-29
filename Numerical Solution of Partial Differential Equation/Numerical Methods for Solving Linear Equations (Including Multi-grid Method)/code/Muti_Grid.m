N = 64; 
h = 1/N;
u = zeros(N-1,1);
f = (sin(pi*[1:N-1]'*h)+sin(16*pi*[1:N-1]'*h))/2;
% 生成矩阵，用于计算残差;
eye1 = 2*ones(1,N-1);
eye2 =-1*ones(1,N-2);
L = 1/h^2*(diag(eye1)+diag(eye2,1)+diag(eye2,-1));
r = zeros(1,11);
r(1) = norm(f-L*u);%初始残差
for i =1:10
u = Muti_Grid_V(f,u,h);
r0 = f-L*u;
r(i+1) = norm(r0);
end


 
v = zeros(N-1,1);
r1 = zeros(1,101);
r1(1) = norm(f-L*v);%初始残差
for i = 1:100
v = sor(L,f,v,1);
r0 = f-L*v;
r1(i+1) = norm(r0);
end
x1 = 0:100;
x2 = 0:10:100;

figure(1)
plot(x2,r,'*-',x1,r1,'+-');
xlabel('迭代次数');
ylabel('$||r_j||_{2}$','Interpreter','latex');
legend('Multi-Grid Method','Gauss-Seidel Method');
title('收敛曲线');

figure(2)
z = zeros(N-1,1);
r_1 = f-L*z;
z = sor(L,f,z,10);
r_2 = f-L*z;
z = sor(L,f,z,90);
r_3 = f-L*z;
r_4 = f-L*u;
x = 0:h:1;
plot(x,[0,r_1',0],'+-',x,[0,r_2',0],'+-',x,[0,r_3',0],'+-',x,[0,r_4',0],'*-');
xlabel('$x_j$','Interpreter','latex');
ylabel('残差r');
legend('0 iterations(G-S)','10 iterations(G-S)','100 iterations(G-S)','100 iterations(Multi-Grid)');
title('残差随迭代次数的变化');




% figure(3)
% u1 = zeros(N-1,1);
% error1 = u1-u;
% u1 = sor(L,f,u1,10);
% error2 = u1-u;
% u1 = sor(L,f,u1,90);
% error3 = u1-u;
% plot(x,[0,error1',0],x,[0,error2',0],x,[0,error3',0]);



function u = sor(L,f,u0,n)
% 进行SOR迭代n次
% w:松弛因子
w = 1;
% u0:迭代初值
D=diag(diag(L)); 
L1=-tril(L,-1); 
U1=-triu(L,1);
u = u0;
for i = 1:n
    u =(D-w*L1)\(((1-w)*D+w*U1)*u+w*f);
end
end