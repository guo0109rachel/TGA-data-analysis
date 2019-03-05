clear all

N=20;
M=15;
dx=1/N;
dy=1/M;

%Mesh
x=linspace(0,1,N+1);
y=linspace(0,1,M+1);

%iteration criteria
e=1e-6;

%elemental matrix
a=(2/dx^2+2/dy^2)*ones(1,N-1);
b_1=-1/dx^2*ones(1,N-2);
b_2=-1/dy^2*ones(1,N-1);

%Matrix construction
MN=(N-1)*(M-1);
A=sparse(MN,MN);
for k=1:M-1
    mm=(N-1)*(k-1)+1;
    nn=(N-1)*k;
    A=A+sparse(mm:nn,mm:nn,a,MN,MN);
    A=A+sparse(mm:(nn-1),(mm+1):nn,b_1,MN,MN)+sparse((mm+1):nn,mm:(nn-1),b_1,MN,MN);
    if k==1
        A=A+sparse((mm+N-1):(nn+N-1),mm:nn,b_2,MN,MN);
    else if k==M-1
            A=A+sparse((mm-N+1):(nn-N+1),mm:nn,b_2,MN,MN);
        else
            A=A+sparse((mm+N-1):(nn+N-1),mm:nn,b_2,MN,MN);
            A=A+sparse((mm-N+1):(nn-N+1),mm:nn,b_2,MN,MN);
        end
    end
end

%f(x,y)=1+x+sin(y);
xf=x(:,2:N);
yf=y(:,2:M);
f=zeros(N-1,M-1);
for a=1:N-1
    for b=1:M-1
        f(a,b)=1+xf(a)+sin(yf(b));
    end
end
f=reshape(f,(N-1)*(M-1),1);

L=tril(A,-1);
D=diag(diag(A,0));
U=triu(A,1);

J=D;%for Jacobi method
GS=L+D;%for gauss-seidal method

%exact solution
ue=A\f;

%Initial condition
v0=zeros(N-1,M-1);
v0=reshape(v0,(N-1)*(M-1),1);

%iteration calculation
nJ=0;
vj=v0;
while max(abs(vj-ue))>e
    vj=vj+D\(f-A*vj);
    nJ=nJ+1;
end

nGS=0;
vgs=v0;
while max(abs(vgs-ue))>e
    %vgs=GS\(f-U*vgs);
    vgs=vgs+GS\(f-A*vgs);
    nGS=nGS+1;
end

V=[ue';vj';vgs'];





