%GS迭代法程序, mseid.m
function [x,k,err,time]=mseidel(A,b,x,tol,max_it)
if nargin<5
    max_it=1000; 
end
if nargin<4
    tol=1.e-5; 
end
if nargin<3
    x=zeros(size(b)); 
end
tic;  
bnrm2=norm(b);
r=b-A*x;  %计算初始残差r0=(b-Ax)
err=norm(r)/bnrm2;
if(err<tol)
    return; 
end
DL=tril(A); 
U=-triu(A,1);
for k=1:max_it    % 迭代开始
    x=DL\(U*x+b);
    r = b-A*x;   %计算残差r=(b-Ax)
    err = norm(r)/bnrm2;
    if (err<=tol)
        break; 
    end
end
time=toc;