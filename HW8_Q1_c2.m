clear all
%set x
N=200;
dx=2/N;
x=linspace(-1,1,N+1);
%set t
M=50;
tf=0.5;
dt=tf/M;
%Problem b1
u1=0.8;
u2=0;
%Set initial Condition
v0=(u1+u2)/2+(u2-u1)/2*tanh(20*x);
v=zeros(M+1,N+1);
v(1,:)=v0;

for n=1:M
    v(n+1,1)=u1;
    v(n+1,N+1)=u2;
    for j=2:N
        v(n+1,j)=1/2*(v(n,j-1)+v(n,j+1))-dt/dx/2*(v(n,j+1)*(1-v(n,j+1))-v(n,j-1)*(1-v(n,j-1)));
    end
end
plot(x,v(M+1,:))