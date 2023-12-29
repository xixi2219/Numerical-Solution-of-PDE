function u = Muti_Grid_V(f,u,h)
% V型多重网格方法,一直把网格变粗到最粗的情况(n=2)再返回
% L:系数矩阵L;
% f:方程右端向量 Lu = f，N-1维;
% u0:迭代初值u_old,N-1维，不包含边界
% h:步长

N = 1/h;
% 生成迭代矩阵L(一维方程情况)
eye1 = 2*ones(1,N-1);
eye2 =-1*ones(1,N-2);
L = 1/h^2*(diag(eye1)+diag(eye2,1)+diag(eye2,-1));
% 在细网格上迭代n1次（光滑迭代）
n1 = 4;
u = sor(L,f,u,n1);
% 计算残差
r = f-L*u;

% 将残差（余量）限制到细网格
rhs = restriction(r);
eps = zeros(size(rhs));

% stop recursion at smallest grid size, otherwise continue recursion
if length(eps) == 1
    eps = 1/2*rhs;
else        
    eps = Muti_Grid_V(rhs,eps,2*h);        
end

 % Prolongation and Correction
 u = u + prolongation(eps);
 
 % Post-Smoothing
 u = sor(L,f,u,n1); 

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
N = (length(r)+1)/2;
res = zeros(N-1,1);
for j = 1:N-1
     res(j) = (r(2*j-1)+2*r(2*j)+r(2*j+1))/4;
end
end

function res = prolongation(eps)
    N = (length(eps)+1)*2;
    res = zeros(N-1,1);
    for j = 2:2:N-2
      res(j) = eps(j/2);  
    end
    for j = 3:2:N-3    
      res(j) = (eps((j+1)/2)+eps((j-1)/2))/2;
    end
      res(1) = eps(1)/2;
      res(N-1) = eps(length(eps))/2;
end