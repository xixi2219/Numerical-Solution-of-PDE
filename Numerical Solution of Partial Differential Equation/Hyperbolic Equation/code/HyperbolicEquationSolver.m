function [u_T,T_max] = HyperbolicEquationSolver(u0,h,meshrate,T,i)
    % Here we consider the initial value problem with perdolical boundary ...
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
% Peridical Condition
a = 1;% the coefficient a, here a>=0
N = length(u0);
tau = meshrate*h/a;% the time step

% 构造更新矩阵
i = [1:N,2:N]; 
j =[1:N,1:N-1];
v = [(1-meshrate)*ones(1,N),meshrate*ones(1,N-1)];
A = sparse(i,j,v,N,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_boundary = [meshrate*u_n(N-1),zeros(1,N-1)]';
    u_next = A*u_n+ u_boundary;
    u_n = u_next;
    T_max = i*tau;
end
 u_T = u_next';
end

% i=2,Lax-Wendroff格式
function [u_T,T_max] = LaxWendroffMethod(u0,h,meshrate,T)
% Peridical Condition
a = 1;% the coefficient a
N = length(u0);
tau = meshrate*h/a;% the time step

% 构造更新矩阵
i = [1:N-1,1:N,2:N]; 
j = [2:N,1:N,1:N-1];
v = [-1/2*meshrate*(1-meshrate)*ones(1,N-1),(1-meshrate.^2)*ones(1,N),1/2*meshrate*(1+meshrate)*ones(1,N-1)];
A = sparse(i,j,v,N,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_boundary = [1/2*meshrate*(1+meshrate)*u_n(N-1),zeros(1,N-2),-1/2*meshrate*(1-meshrate)*u_n(2)]';
    u_next = A*u_n+ u_boundary;
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
i = [1:N,2:N,3:N]; 
j = [1:N,1:N-1,1:N-2];
v = [1/2*(1-meshrate)*(2-meshrate)*ones(1,N),meshrate*(2-meshrate)*ones(1,N-1),-1/2*meshrate*(1-meshrate)*ones(1,N-2)];
A = sparse(i,j,v,N,N);

% 迭代计算数值解
TimeStepNum = T/tau;% 迭代次数
u_next = zeros(length(u0),1); 
for i=1:TimeStepNum
    if i==1
        u_n = u0';
    end
    u_boundary = [meshrate*(2-meshrate)*u_n(N-1)-1/2*meshrate*(1-meshrate)*u_n(N-2),-1/2*meshrate*(1-meshrate)*u_n(N-1),zeros(1,N-2)]';
    u_next = A*u_n+ u_boundary;
    u_n = u_next;
    T_max = i*tau;
end
 u_T = u_next';
end