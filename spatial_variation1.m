x=linspace(0,1,21);
x=x';

ui=zeros(19,1);
ui(10)=1;
a=spatial_var1(1,ui,5,5000,20);

%{
ui=sin(pi/2*x);
a=spatial_var1(1,ui,5,5000,21);
%}

function [g]=spatial_var1(K,u0,tn,n,N)  %u0 is a vector in [0,1], N is number of grid points in use

u=zeros(N-1,n);  %store all u
u(:,1)=u0;  %assume boundary points are fictitous and used to determine first and last entries in the matrix
int=zeros(N-1,1);   %store integrals
dx=1/N;
x=dx:dx:(1-dx);
t=linspace(0,tn,n);
dt=tn/(n-1);
s=dt/K/(dx)^2;   

A=full(gallery('tridiag',N-1,s,1-2*s,s));  %create matrix
A(1,1)=1-s;
A(N-1,N-1)=1-s;  %neumann boundary condition

for i=1:(n-1)
   if i==1
       int=zeros(N-1,1);
   else
       int=int+dt*(u(:,i)+u(:,i-1))/2;
   end
   
   u(:,i+1)=A*u(:,i)+dt/K*u(:,i)-dt/K*u(:,i).*u(:,i)-dt/K*u(:,i).*int;
   bar(x,u(:,i+1));
   pause( 0.05 );
end
g=int;
end