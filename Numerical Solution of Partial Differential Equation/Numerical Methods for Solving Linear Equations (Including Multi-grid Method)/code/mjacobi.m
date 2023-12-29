%Jacobi����������  mjacobi.m
 function [x,k,err,time]=mjacobi(A,b,x,tol,max_it)
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
 bnrm2 = norm(b);
 r=b-A*x;  %�����ʼ�в�r0=(b-Ax)
 err=norm(r)/bnrm2;
 if (err<tol)
     return; 
 end
 D=diag(diag(A));
 for k=1:max_it   % ������ʼ
     x=D\((D-A)*x+b);
     r=b-A*x;   %����в�r=(b-Ax)
     err=norm(r)/bnrm2;
     if(err<=tol)
         break; 
     end
 end
 time=toc;