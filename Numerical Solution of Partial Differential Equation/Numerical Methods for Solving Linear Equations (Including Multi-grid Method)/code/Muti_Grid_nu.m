N = 64; 
h = 1/N;
u = zeros(N-1,1);
f = (sin(pi*[1:N-1]'*h)+sin(16*pi*[1:N-1]'*h))/2;
% 生成矩阵，用于计算残差;
eye1 = 2*ones(1,N-1);
eye2 =-1*ones(1,N-2);
L = 1/h^2*(diag(eye1)+diag(eye2,1)+diag(eye2,-1));
i = 0;%指标
r = f-L*u;
tic;
while norm(r)>1e-10
    u = Muti_Grid_V(f,u,h);
    r = f-L*u;
    i = i+1;
end
time = toc;
r0 = norm(r);