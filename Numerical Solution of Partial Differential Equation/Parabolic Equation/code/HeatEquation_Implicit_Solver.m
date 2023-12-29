function [Error1,Error2,varargout]=HeatEquation_Implicit_Solver(N,mu)
%% 1-D problem Implicit Method for smooth initial value

h = 1/N; % Space Steps
x = linspace(0,1,N+1); % grid point
T0 = 1;

Delta_t = h^2 * mu;
NumT = T0/ Delta_t; % the Num of T Steps
U_0 = sin(pi*x(2:end-1)); % Initial Condition

% Generate the Coff Matrix 
% mu = Delta_t/h^2;
Main_Diag = (1+2*mu)*ones(1,N-1);
Diag_U1 = -mu * ones(1,N-1);
Diag_D1 = -mu * ones(1,N-1);
% A = diag(Main_Diag,0) + diag(Diag_D1,-1) + diag(Diag_U1,1);
A = spdiags(Main_Diag',0,N-1,N-1);
A = spdiags(Diag_U1',1, A);
A = spdiags(Diag_D1',-1,A);


for n = 0 : 1 : NumT-1
    if n == 0
       U_n = U_0'; 
    end
 
    U_next = A\U_n;

    U_n = U_next;
    M = n;
end

% Exact Value
T1 = M*Delta_t;
Exact_U = exp(-pi^2 * T1)*sin(pi*x(2:end-1));
Error1 = max(abs(Exact_U'-U_n));
Error2 =sqrt( 1/(N) * sum(abs(Exact_U'-U_n).^2));

if nargout == 3
    varargout{1} = U_n;
end

end