function u = Muti_Grid_V_2D(f,u,h)
% 2D情况下的V型多重网格方法,一直把网格变粗到最粗的情况(n=2)再返回
% L:系数矩阵L;
% f:方程右端向量 Lu = f，（N-1）*（N-1）维;
% u:迭代初值u_old,（N-1）*（N-1）维，不包含边界
% h:步长

N = 1/h;
% 生成迭代矩阵L(二维方程情况)
n = N-1;%矩阵块维数,最终所求得的离散矩阵为(N-1)^2阶
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
L = (1/h^2)*kron(D,I)+kron(I,D);%Kronecker运算求得离散矩阵


% 在细网格上迭代n1次（光滑迭代）
n1 = 1;
% 对u做重新排序
u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
u1 = sor(L,f,u1,n1);
u(2:N,2:N) = reshape(u1,[N-1,N-1]);
u(2:N,2:N) = u(2:N,2:N)';

% 计算残差
r1 = f-L*u1;
r = u;
r(2:N,2:N) = reshape(r1,[N-1,N-1]);
r(2:N,2:N) = r(2:N,2:N)';

% 将残差（余量）限制到细网格
rhs = restriction(r);
eps = zeros(size(rhs));

% stop recursion at smallest grid size, otherwise continue recursion
if length(eps)-1 == 2
    eps(2,2) = 1/4*rhs(2,2);
else
    rhs = rhs(2:end-1,2:end-1);
    rhs = rhs';
    rhs = reshape(rhs,[length(rhs)*length(rhs),1]); 
    eps = Muti_Grid_V_2D(rhs,eps,2*h);        
end

 % Prolongation and Correction
 u = u + prolongation(eps);
 
 % Post-Smoothing
 u1 = reshape(u(2:N,2:N)',[(N-1)*(N-1),1]);
 u1 = sor(L,f,u1,n1); 
 u(2:N,2:N) = reshape(u1,[N-1,N-1]);
 u(2:N,2:N) = u(2:N,2:N)';
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

function res = restriction(r)
% 九点加权平均
N = (length(r)-1)/2;
res = zeros(N+1,N+1);
for i = 2:N
    for j = 2:N
     res(i,j) = (4*r(2*i-1,2*j-1)+2*r(2*i-2,2*j-1)+2*r(2*i,2*j-1)+2*r(2*i-1,2*j-2)+2*r(2*i-1,2*j)+r(2*i-2,2*j-2)+r(2*i-2,2*j)+r(2*i,2*j-2)+r(2*i,2*j))/16;
    end
end
end

function res = prolongation(eps)
% 九点延拓
    N = (length(eps)-1)*2;
    res = zeros(N+1,N+1);
    for i = 1:2:N+1
        for j = 1:2:N+1
            res(i,j) = eps((i+1)/2,(j+1)/2);
        end
        for j = 2:2:N
            res(i,j) = (eps((i+1)/2,j/2)+eps((i+1)/2,j/2+1))/2;
        end
    end
    
    
    for i  = 2:2:N
        for j = 1:2:N+1
            res(i,j) = (eps(i/2,(j+1)/2)+eps(i/2+1,(j+1)/2))/2;
        end
        for j = 2:2:N
            res(i,j) = (res(i-1,j-1)+res(i-1,j+1)+res(i+1,j+1)+res(i+1,j-1))/4;
        end
    end 
end