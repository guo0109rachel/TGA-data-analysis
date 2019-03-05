%Calculation of fibrous medium radiative and conductive heat transfer
%Assumed wavelength: ?=4.0 um
wl=4e-6;
%Medium density: 20 kg/m3
rou=20;
%Thickness: E=10 cm
E=0.1;

%constants for thermal conductivity
a=0.2572;
c=0.0527*rou^0.91;
b=0.0013*c;

%Temperature: T_0=400 K,T_E=300 K
T0=400;
TE=300;

psi=@(x)a/1.81*x^1.81+b/2*x^2+c*x;
Tls0=psi(T0);
TlsE=psi(TE);

%Spectral coefficient
%absorption 500 m^-1
sa=500;
%scattering 1000 m^-1
ss=1000;
%extinction 1500 m^-1
se=sa+ss;

%Mesh
N=100;%mesh for position
y=linspace(0,E,N+1);
y=y(:,2:N);
h=E/N;

m=12;%mesh for angular
u=linspace(1,0,m/2+2);
u=u(:,2:(m/2+1));

%Boundary condition
L0=BBMI(wl,T0)*ones(1,m/2);
LE=BBMI(wl,TE)*ones(1,m/2);

%Initial temperature profile guess
TI=linspace(T0,TE,N+1);
T=TI(:,2:N);

%Newton-Cotes integration quadrature weights
QW5=NCQWO5();
C=[QW5,QW5];
%QW14=NCQW14();
%C=[QW14(14),QW14(13),QW14(12),QW14(11),QW14(10),QW14(9),QW14(:,2:7)];

%Monochromatic coefficients for scattering calculation
Cm=MCco();

%Iteration
P=3;
T_err=zeros(P,N-1);
Ts=zeros(P,N-1);
for p=1:P

%matrix A_lambda
%Alb=zeros(m,m);
Alb=Cm-se*diag(ones(1,m));
ur=[1./u,1./flip(-u)];
Alb=diag(ur)*Alb;
%for j=1:m/2
%   Alb(j,:)=1/u(j)*sP*ones(1,m);
%   Alb(j,j)=Alb(j,j)-se/u(j);
%   Alb(m/2+j,:)=-1/u(j)*sP*ones(1,m);
%   Alb(j+m/2,j+m/2)=Alb(j+m/2,j+m/2)+se/u(j);
%end

%G matrix
%sub G matrices
G2=-2*h*Alb;
G4=diag([ones(1,m/2),4*ones(1,m/2)]);
G5=diag([zeros(1,m/2),-ones(1,m/2)]);
G6=diag([ones(1,m/2),zeros(1,m/2)]);
G7=diag([-4*ones(1,m/2),-ones(1,m/2)]);
G8=diag([zeros(1,m/2),-3*ones(1,m/2)]);
G9=diag([3*ones(1,m/2),zeros(1,m/2)]);
G1=G2+G8;
G3=G2+G9;

%construct matrix G
M=(N-1)*m;
G=sparse(M,M);
for j=1:m
    G=G+sparse(j*ones(1,m),1:m,G1(j,:),M,M);
    G=G+sparse(j*ones(1,m),m+(1:m),G4(j,:),M,M);
    G=G+sparse(j*ones(1,m),2*m+(1:m),G5(j,:),M,M);
    G=G+sparse(M-m+j*ones(1,m),M-m+(1:m),G3(j,:),M,M);
    G=G+sparse(M-m+j*ones(1,m),M-2*m+(1:m),G7(j,:),M,M);
    G=G+sparse(M-m+j*ones(1,m),M-3*m+(1:m),G6(j,:),M,M);
end

for i=2:N-2
    G=G+sparse((i-1)*m+(1:m),(i-2)*m+(1:m),-ones(1,m),M,M);
    G=G+sparse((i-1)*m+(1:m),i*m+(1:m),ones(1,m),M,M);
    for j=1:m
        G=G+sparse((i-1)*m+j*ones(1,m),(i-1)*m+(1:m),G2(j,:),M,M);
    end
end

%Matrix F
F=zeros(N-1,m);
F(1,:)=2*h*[sa*BBMI(wl,T(1))./u,-sa*BBMI(wl,T(1))./u]+[L0,zeros(1,m/2)];
F(N-1,:)=2*h*[sa*BBMI(wl,T(N-1))./u,-sa*BBMI(wl,T(N-1))./u]+[zeros(1,m/2),LE];
for i=2:N-2
    F(i,:)=2*h*[sa*BBMI(wl,T(i))./u,-sa*BBMI(wl,T(i))./u];
end
F=reshape(F',M,1);

%solve matrix K
K=G\F;

%radiative heat flux Qr
Qr=zeros(1,N-1);
for i=1:N-1
    cut=((i-1)*m+1):i*m;
    Qr(i)=2*pi()*sum(transpose(K(cut,:)).*[u,-flip(u)].*C);
end

%radiative source term Sr
Sr=zeros(N-1,1);
Sr(1)=-1/12/h*(-25*Qr(1)+48*Qr(2)-36*Qr(3)+16*Qr(4)-3*Qr(5));
Sr(2)=-1/12/h*(-3*Qr(1)-10*Qr(2)+18*Qr(3)-6*Qr(4)+Qr(5));
Sr(N-2)=-1/12/h*(-Qr(N-5)+6*Qr(N-4)-18*Qr(N-3)+10*Qr(N-2)+3*Qr(N-1));
Sr(N-1)=-1/12/h*(3*Qr(N-5)-16*Qr(N-4)+36*Qr(N-3)-48*Qr(N-2)+25*Qr(N-1));
for i=3:N-3
    Sr(i)=-1/12/h*(Qr(i-2)-8*Qr(i-1)+8*Qr(i+1)-Qr(i+2));
end

%step 3.1
D=diag(-2*ones(1,N-1),0);
D=D+diag(ones(1,N-2),1)+diag(ones(1,N-2),-1);
Sr=-h^2*Sr-[Tls0;zeros(N-3,1);TlsE];
Tls=D\Sr;

%function 36
Tn=zeros(1,N-1);
for i=1:N-1
    FKT=@(x)a/1.81*x^1.81+b/2*x^2+c*x-Tls(i);
    Tn(i)=fsolve(FKT,350);
end

T_err(p,:)=(Tn-T)./Tn;
%T_err(p,:)=real(Tn-T);
T=Tn;
Ts(p,:)=Tn;
%T=real(Tn);

end