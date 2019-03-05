function f = fun( t, x )
km=0.01725;
Ea=2660;
R=1;
Tm=300;
Ca0=2;
Q=-41700;
CT=80;
V=1200;
H=-10000;
k=km*exp(-Ea/R.*(1./x(2)-1/Tm));
f=zeros(2,1);
f(1)=k*Ca0*(1-x(1)).^2;
f(2)=Q/CT/V-H*Ca0/CT*f(1);
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

