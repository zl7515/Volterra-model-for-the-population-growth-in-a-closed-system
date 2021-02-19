ui=zeros(19,19);
ui(10,10)=1;
spatial_var2(1,ui,1,5000,20);

function [g]=spatial_var2(K,u0,tn,n,N)      %u0 is an N-1*N-1 matrix, N is number of mesh in each direction

u=zeros((N-1)*(N-1),n);
k=1;
mesh=1/N;
x=mesh:mesh:(1-mesh);
dt=tn/(n-1);
int=zeros((N-1)*(N-1),1);
s=dt/K/mesh^2;

%change matrix to vector
for i=1:N-1
    for j=1:N-1
        u(k,1)=u0(i,j);
        k=k+1;
    end
end

A=makeA(N);
A=-A;

%neumann boundary condition
for i=1:(N-1)*(N-1)
    row_sum=sum(A(i,:));
    if row_sum==-2
        A(i,i)=-2;
    elseif row_sum==-1
        A(i,i)=-3;
    end
end
        
for i=1:(n-1)
   if i==1
       int=zeros((N-1)*(N-1),1);
   else
       int=int+dt*(u(:,i)+u(:,i-1))/2;
   end
   
   u(:,i+1)=s*A*u(:,i)+(1+dt/K)*u(:,i)-dt/K*u(:,i).*u(:,i)-dt/K*u(:,i).*int;
   u_temp=u(:,i+1);
   
   %restore to matrix form
   u_mat=zeros(N-1,N-1);
   ind=1;
   for j=1:N-1
       for k=1:N-1
           u_mat(j,k)=u_temp(ind);
           ind=ind+1;
       end
   end
   
   %visualization
   
   surf(x,x,u_mat);
   pause( 0.05 );
end

end