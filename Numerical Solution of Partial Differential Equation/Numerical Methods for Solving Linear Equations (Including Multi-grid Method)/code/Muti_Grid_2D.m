N = 64; 
h = 1/N;
x = 0:h:1; y = 0:h:1;
[X,Y] = meshgrid(x,y);%格点坐标
u = zeros(N+1,N+1);
f =  Trigonometric_Func(x , y);
f = f(2:end-1,2:end-1);
f = f';
f = reshape(f,[(N-1)*(N-1),1]); % f的值

% 生成迭代矩阵L(二维方程情况)
n = N-1;%矩阵块维数,最终所求得的离散矩阵为(N-1)^2阶
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
L = (1/h^2)*kron(D,I)+kron(I,D);%Kronecker运算求得离散矩阵

r = zeros(11,1);
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
r(1) = max(abs(L*u1-f)); 
for i = 1:10
   u = Muti_Grid_V_2D(f,u,h);
   u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
   r0 = L*u1-f;
   r(i+1) = max(abs(r0));  
end

v = zeros(N+1,N+1);
v1 = reshape(v(2:N,2:N)',[(N-1)*(N-1),1]);
r1 = zeros(101,1);
r1(1) = max(abs(L*v1-f)); 
for i = 1:100 
    v1 = sor(L,f,v1,1);
    r0 = L*v1-f;
    r1(i+1) = max(abs(r0));     
end
figure(1)
x1 = 0:100;
x2 = 0:10:100;
plot(x2,r,'*-',x1,r1,'+-');
xlabel('迭代次数');
ylabel('$||r_j||_{\infty}$','Interpreter','latex');
legend('Multi-Grid Method','Gauss-Seidel Method');
title('收敛曲线');


figure(2)
% 一次多重网格
subplot(1,3,1)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
u = Muti_Grid_V_2D(f,u,h);
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('1 iteration');

% 十次多重网格
subplot(1,3,2)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
for i =1:10
u = Muti_Grid_V_2D(f,u,h);
end
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('10 iterations');

% 100次多重网格
subplot(1,3,3)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
for i =1:100
u = Muti_Grid_V_2D(f,u,h);
end
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('100 iterations');

figure(3)
% 初始
subplot(1,3,1)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('0 iteration');

% 10次G-S
subplot(1,3,2)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
u1 = sor(L,f,u1,10);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('10 iterations');

% 1000次G-S
subplot(1,3,3)
u = zeros(N+1,N+1);
r_1= zeros(N+1,N+1);
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
u1 = sor(L,f,u1,100);
r0 = L*u1-f;
r_1(2:N,2:N) = reshape(r0,[N-1,N-1]);
r_1(2:N,2:N) = r_1(2:N,2:N)';
mesh(X,Y,r_1);
title('100 iterations');
function [f] =  Trigonometric_Func(x,y)
% -f在网格点上的值
    [X,Y] = meshgrid(x,y);
    f = sin(pi*X).*cos(pi*Y)+sin(16*pi*X).*cos(16*pi*Y);
end

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