function [Error1,Error2,varargout]=HeatEquation_CN_Solver(N,mu)
%% 1-D problem CN Method for smooth initial value

h = 1/N; % Space Steps
x = linspace(0,1,N+1); % grid point
T0 = 1;

Delta_t = h^2*mu;
NumT = T0/ Delta_t; % the Num of T Steps
U_0 = sin(pi*x(2:end-1)); % Initial Condition

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
    M = n;
end

% Exact Value
T1 = M*Delta_t;
Exact_U = exp(-pi^2*T1)*sin(pi*x(2:end-1));
Error1 = max(abs(Exact_U'-U_n));
Error2 =sqrt( 1/(N) * sum(abs(Exact_U'-U_n).^2));

if nargout == 3
    varargout{1} = U_n;
end
end