function [U,e] = PossionEquation1(N,F,M,u)
% 在Dirichlet边界条件下求解方程: -\Delta u= f  
% U 数值解 
% X x轴坐标,Y y轴坐标;
% e 误差的最大范数
% F -f在网格点上的值
% M 边界值
% u 真解在网格点上的值

if (size(F,1) ~= N+1) 
   return
end

h= 1/N; %分割段数
x = 0:h:1; y = 0:h:1;
[X,Y] = meshgrid(x,y);%格点坐标

%% laplace离散矩阵的构造
n = N-1;%矩阵块维数,最终所求得的离散矩阵为(N-1)^2阶
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
L = kron(D,I)+kron(I,D);%Kronecker运算求得离散矩阵

%% Dirchlet边界条件的构造 1.uD=sin(pi*x)cos(pi*y) 2.uD=x^2+y^2
M(2:end-1,2:end-1) = 0;
M(2:end-1,2) = M(2:end-1,2)+M(2:end-1,1);
M(2:end-1,end-1) = M(2:end-1,end-1)+M(2:end-1,end);
M(2,2:end-1)= M(2,2:end-1) + M(1,2:end-1);
M(end-1,2:end-1) = M(end-1,2:end-1) + M(end,2:end-1);
M = M(2:end-1,2:end-1);

M = M';
M = reshape(M,[(N-1)*(N-1),1]); % 边界条件

%% 函数-f的计算1.-f=-2*pi^2*sin(pi*x)cos(pi*y) 2. -f = -(x^2+y^2)exp(x+y)
F = F(2:end-1,2:end-1);
F = F';
F = reshape(F,[(N-1)*(N-1),1]); % f的值


%%  求解线性方程组
U = (L)\(h^2 * F-M);

U = reshape(U,[N-1,N-1]);
U = U';

%% 真解
if nargin ==4 & nargout==2

   u = u(2:end-1,2:end-1);
    
   e=max(max(abs(U-u)));
    
end

end