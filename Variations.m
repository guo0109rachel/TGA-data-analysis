clear all
%Set parameters
N=20;
h=1/N;
dt=1/20;
Tstep=20;
r=dt/h^2;

N1=N+1;
N2=N1^2;

%set x and y on the grid
x=linspace(0,1,N1);
y=linspace(0,1,N1);

%Set initial conditions
%v=1+x+x.^2+y+y.^2;
v=zeros(N1);
for k=1:N1
    for j=1:N1
        v(j,k)=1+h*(j-1)+(h*(j-1))^2+h*(k-1)+(h*(k-1))^2;
    end
end

%Build Martix A
A=eye(N2); %BCs at y=0 & y=1

for j=(N1+1):(N2-N1)
    if rem(j-1,N1)==0 %BCs at x=0
    else if rem(j,N1)==0 %BCs at x=1
    else A(j,j)=1+2*r;
         A(j,j-1)=-r/2;
         A(j,j+1)=-r/2;
         A(j,j-N1)=-r/2;
         A(j,j+N1)=-r/2;
        end
    end
end
% LU factorization
[L,U,P]=lu(A);

%Set Cn
Cn=zeros(N2,1);
for n=1:20
    %Construct Cn
    for i=1:N1 %BC y=0 & y=1
        Cn(i)=1+h*(i-1)+(h*(i-1))^2+n*dt+(n*dt)^2;
        Cn(i+N2-N1)=3+h*(i-1)+(h*(i-1))^2+n*dt+(n*dt)^2;
    end
    for i=(N1+1):(N2-N1)
        if rem(i-1,N1)==0 %BC x=0
            k=floor(i/N1);
            Cn(i)=1+k*h+(k*h)^2+n*dt+(n*dt)^2;
        else if rem(i,N1)==0 %BC x=1
                k=floor(i/N1)-1;
                Cn(i)=3+k*h+(k*h)^2+n*dt+(n*dt)^2;
            else j=rem(i-1,N1);
                k=floor(i/N1);
                Cn(i)=(1-2*r)*v(j+1,k+1)+r/2*(v(j,k+1)+v(j+2,k+1)+v(j+1,k)+v(j+1,k+2))+1/2*(2*dt*n-3+2*dt*(n-1)-3)*dt;
            end
        end
    end
    %solve the function
    vn=L\(P*Cn);
    vn=U\vn;
    %solution for n+1
    v=reshape(vn,N1,N1);
end
%exact solution
%u=1+x+x.^2+y+y.^2+t+t.^2;
u=zeros(N1);
for k=1:N1
    for j=1:N1
        u(j,k)=3+x(j)+x(j)^2+y(k)+y(k)^2;
    end
end
uerror=u-v;
surf(x,y,uerror);
error= max(abs(uerror(:)));
