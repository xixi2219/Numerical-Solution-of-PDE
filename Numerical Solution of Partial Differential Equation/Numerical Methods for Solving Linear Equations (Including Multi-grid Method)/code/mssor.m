%SSOR迭代法程序, mssor.m
function [x,k,err,time]=mssor(A,b,w,x,tol,max_it)
if nargin<6
    max_it=1000; 
end
if nargin<5
    tol=1.e-6; 
end
if nargin<4
    x=zeros(size(b)); 
end
tic; 
bnrm2 = norm(b);
r=b-A*x;  %计算初始残差r0=(b-Ax)
err=norm(r)/bnrm2;
if(err<tol)
    return; 
end
D=diag(diag(A));
L=-tril(A,-1); 
U=-triu(A,1);
for k=1:max_it    % 迭代开始
    x=(D-w*L)\(((1-w)*D+w*U)*x+w*b);
    x=(D-w*U)\(((1-w)*D+w*L)*x+w*b);
    r=b-A*x;     %计算残差r=(b-Ax)
    err=norm(r)/bnrm2;
    if(err <= tol)
        break; 
    end
end
time=toc;