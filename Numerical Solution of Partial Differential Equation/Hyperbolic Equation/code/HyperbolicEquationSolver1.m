function [u_T,T_max] = HyperbolicEquationSolver1(u0,h,meshrate,T,i)
    % Here we consider the initial value problem with left boundary
    %(此时,第1个点个第N个点即为边界点)
    % u_t + a u_x = 0 . Here a = 1.
    % u0: the initial value.
    % h: the space step.
    % meshrate: the mesh rate defined by a*tau/h.(here tau means the time step)
    % T:t_max,the max time we computate.
    % u_T: 所求数值解
    % T_max: 最终计算到的时间步,用于计算真解
    if i==1
       [u_T,T_max] = UpwindMethod(u0,h,meshrate,T);
    elseif i==2
       [u_T,T_max] = LaxWendroffMethod(u0,h,meshrate,T);
    elseif i==3
       [u_T,T_max] = BeamWarmingMethod(u0,h,meshrate,T);
    else
        error("没有对应的数值方法！")
    end
    
end
% i=1，迎风格式  
function [u_T,T_max] = UpwindMethod(u0,h,meshrate,T)
% zero boundary
a = 1;% the coefficient a, here a>=0
N = length(u0);
tau = meshrate*h/a;% the time step

% 构造更新矩阵
i = [1:N-1,1:N-1]; 
j =[2:N,1:N-1];
v = [(1-meshrate)*ones(1,N-1),meshrate*ones(1,N-1)];
A = sparse(i,j,v,N-1,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_next(2:end) = A*u_n;
    u_next(1) =(1-meshrate)*u_n(1)+meshrate*u_n(1);% 算例使用 zero boundary u0(1) = 0; 零阶外推公式处理边界
    u_n = u_next;
    T_max = i*tau;
end
 u_T = u_next';
end

% i=2,Lax-Wendroff格式
function [u_T,T_max] = LaxWendroffMethod(u0,h,meshrate,T)
% zero boundary
a = 1;% the coefficient a
N = length(u0);
tau = meshrate*h/a;% the time step

% 构造更新矩阵
i = [1:N-2,1:N-2,1:N-2]; 
j = [3:N,2:N-1,1:N-2];
v = [-1/2*meshrate*(1-meshrate)*ones(1,N-2),(1-meshrate.^2)*ones(1,N-2),1/2*meshrate*(1+meshrate)*ones(1,N-2)];
A = sparse(i,j,v,N-2,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_next(2:N-1) = A*u_n;
    u_next(1) = -1/2*meshrate*(1-meshrate)*u_n(1)+(1-meshrate.^2)*u_n(1)+1/2*meshrate*(1+meshrate)*u_n(1);% 算例里使用zero boundary u0(1)=0，但在这里是用零阶外推公式处理边界,也可以输入其它边界
    u_next(N) = -1/2*meshrate*(1-meshrate)*u_n(N)+(1-meshrate.^2)*u_n(N)+1/2*meshrate*(1+meshrate)*u_n(N-1);% 零阶外推公式
    u_n = u_next;
    T_max = i*tau;
end
 u_T = u_next';
end


function [u_T,T_max] = BeamWarmingMethod(u0,h,meshrate,T)
% Peridical Condition
a = 1;% the coefficient a, here a>=0
N = length(u0);
tau = meshrate*h/a;% the time step

% 构造更新矩阵
i = [1:N-2,1:N-2,1:N-2]; 
j = [3:N,2:N-1,1:N-2];
v = [1/2*(1-meshrate)*(2-meshrate)*ones(1,N-2),meshrate*(2-meshrate)*ones(1,N-2),-1/2*meshrate*(1-meshrate)*ones(1,N-2)];
A = sparse(i,j,v,N-2,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_next(3:N) = A*u_n;
    u_next(1) = 1/2*(1-meshrate)*(2-meshrate)*u_n(1)+meshrate*(2-meshrate)*u_n(1)-1/2*meshrate*(1-meshrate)*u_n(1);% 算例里使用zero boundary u0(1)=0，但在这里是用零阶外推公式处理边界,也可以输入其它边界
    u_next(2) = 1/2*(1-meshrate)*(2-meshrate)*u_n(2)+meshrate*(2-meshrate)*u_n(1)-1/2*meshrate*(1-meshrate)*u_n(1);
    u_n = u_next;
    T_max = i*tau;
end
 u_T = u_next';
end