function varargout=HeatEquation_CN_Solver_1(N,mu,U_0)
%% 1-D problem CN Method for unsmooth initial value

h = 1/N; % Space Steps
x = linspace(0,1,N+1); % grid point
T0 = 1;

Delta_t = h^2*mu;
NumT = T0/ Delta_t; % the Num of T Steps

% Generate Left Matrix and Right Matrix
% left Matrix
% mu = Delta_t/h^2;
Main_Diag = (1+mu) * ones(1,N-1);
Diag_U1 = (-mu/2) * ones(1,N-1);
Diag_D1 = (-mu/2) * ones(1,N-1);
A = spdiags(Main_Diag',0,N-1,N-1);
A = spdiags(Diag_U1',1, A);
A = spdiags(Diag_D1',-1,A);
% Right Matrix
Main_Diag = (1-mu) * ones(1,N-1);
Diag_U1 = (mu/2) * ones(1,N-1);
Diag_D1 = (mu/2) * ones(1,N-1);
B = spdiags(Main_Diag',0,N-1,N-1);
B = spdiags(Diag_U1',1, B);
B = spdiags(Diag_D1',-1,B);

for n = 0:1:NumT-1
   if n == 0
      U_n = U_0'; 
   end
   
    b = B* U_n;
   
    U_next = A \ b;
    U_n = U_next;
end
    varargout{1} = U_n;
end