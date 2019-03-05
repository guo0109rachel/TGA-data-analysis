function uout=pdepa(tlist,xlist,u0,bd,T,k,X,h,ax,bx,cx,dx,theta)
% Simple solver for a 1-dimensional linear parabolic PDE.
% Use a weighted implicit scheme.
% Parameters:
% tlist is the time corresponding to desired results;
% xlist is the space point corresponding to desired results;
% u0 is the initial condition, either a constant or a function handle, corresponding to each of x;
% bd is a CELL of boundary conditions;
% T contains the boundaries of time grid, and k contains the difference
% between two grid points.
% They should be an array in the form T=[T0,T1,T2,...], and k=[k1,k2,...]
% Time between T0 and T1 is discretized with step sizes of k1, and time
% between T1 and T2 with step sizes of k2.
% X and h describe the grid of space, in the format similar to T and k:
% X=[X0,X1,X2,...], and h=[h1,h2,...]
% ax ,bx and cx are parameters in the equation:
% du/dt=a*d2u/dx2+b*du/dx+c*u+d
% they must be independent of time, and should be given in the form of
% scalar or function handle.
% theta is the weight number. It should be in the interval [0,1].
% theta=0: classical explicit scheme; theta=0.5: Crank-Nicolson scheme;
% theta=1: classical implicit scheme.
uout = zeros(length(tlist),length(xlist));
uout1 = zeros(1,length(xlist));
uout2 = uout1;
% Grid of space and time
x=X(1):h(1):X(2);
xbound=length(x); % xbound contains indexes of x interval bounds.
if length(X)>2
for i=2:length(X)
x=[x,X(i)+h(i):h(i):X(i+1)];
xbound(i)=length(x);
end
end
hx=diff(x);
l=length(x);
t=T(1):k(1):T(2);
tbound=length(t); % tbound contains indexes of t interval bounds.
if length(T)>2
for i=2:length(T)
t=[t,T(i)+k(i):k(i):T(i+1)];
tbound(i)=length(t);
end
end
kt=diff(t);
% expend parameters a, b and c:
if isnumeric(ax)
a=ax+zeros(1,l);
else
a=ax(x);
end
if isnumeric(bx)
b=bx+zeros(1,l);
else
b=bx(x);
end
if isnumeric(cx)
c=cx+zeros(1,l);
else
c=cx(x);
end
if isnumeric(dx)
d=dx+zeros(1,l);
else
d=dx(x);
end
% % Stability
% if theta<0.5 && max(a)*max(k)/min(h)^2>0.5/(1-2*theta)
% fprintf('Not stable !')
% return;
% end
n_tint=0;% the number of time intervals calculated, or basic equations constructed.
n_tout=1;
if isnumeric(u0)
u1=u0+zeros(1,l);
else
u1=u0(x);
end
u1=u1';
% In the NO.tstep circle, calculate the value at time t(tstep+1)
for tstep=1:length(t)-1
if tstep == 1 || tstep == tbound(n_tint)
n_tint=n_tint+1;
% Basic equations: A*u(t+1)=B*u(t)+C
% Basic equations, Left side:
% Row m is an equation of the point m+1; column n corresponds to point n.
% The last two rows are boundary conditions.
A=sparse(2:l-1,1:l-2,-1.*(a(2:l-1)./hx(1:l-2)-b(2:l-1)./2).*theta.*k(n_tint)./hx(1:l-2),l,l);
A=A+sparse(2:l-1,2:l-1,(2*a(2:l-1)./hx(1:l-2).^2-c(2:l-1)).*theta.*k(n_tint)+1,l,l);
A=A+sparse(2:l-1,3:l,-1.*(a(2:l-1)./hx(1:l-2)+b(2:l-1)./2).*theta.*k(n_tint)./hx(1:l-2),l,l);
% Basic equations, Right side:
B=sparse(2:l-1,1:l-2,(a(2:l-1)./hx(1:l-2)+b(2:l-1)./2).*(1-theta).*k(n_tint)./hx(1:l-2),l,l);
B=B+sparse(2:l-1,2:l-1,(-2*a(2:l-1)./hx(1:l-2).^2+c(2:l-1)).*(1-theta).*k(n_tint)+1,l,l);
B=B+sparse(2:l-1,3:l,(a(2:l-1)./hx(1:l-2)-b(2:l-1)./2).*(1-theta).*k(n_tint)./hx(1:l-2),l,l);
C=sparse(2:l-1,1,k(n_tint).*d(2:l-1),l,1);
% Default boundary conditions at interval bounds: continuous derivative.
if xbound(1)<l
for i=1:length(xbound)-1
j=xbound(i);
A(j-1,j-2:j+2)=[1,-4,3,0,0]/h(i)-[0,0,-3,4,-1]/h(i+1);
B(j-1,:)=0;
C(j-1,1)=0;
end
end
% Boundary conditions: do not vary with time
% Format of bd: {type, x, value}. bd is a cell!!
% types: 1: boundary function value;
% 2: left derivative value, used in the right boundary of a region;
% 3: right derivative value, used in the left boundary of a region;
% 4: relation between two derivative at a bound(conductive boundary condition);
% 5,6,7: custumer defined.
% type 5,6,7 boundary condition should be a cell, in the format: {type, x, value, formula}
% 'formula' is a function handle, with parameter h and returns a row array of odd length.
% 'formula' can also be a row array.
% For type 5, the last term of 'formula' corresponds to x;
% for type 6, the first term corresponds to x;
% for type 7, the central term corresponds to x.
n=0; % number of boundary conditions added.
for i=1:size(bd,1)
j=find(x==bd{i,2});
if bd{i,1}~=4
if j==1 || j==l
C(j,1) = bd{i,3};
else
C(l-1+n,1)=bd{i,3};
end
end
switch bd{i,1}
case 1 % boundary function value
if j==1 || j==l
A(j,j) = 1;
else
A(l-1+n,j)=1;
end
case 2 % left derivative
A(l-1+n,j-2:j)=[0.5,-2,1.5]./hx(j-1);
case 3 % right derivative
A(l-1+n,j:j+2)=[-1.5,2,-0.5]./hx(j);
case 4 % du/dx(left)-k[du/dx(right)]=0; k is 'value' in bd.
A(j-1,j-2:j+2)=[1,-4,3,0,0]-bd{i,3}*[0,0,-3,4,-1];
B(j-1,:)=0;
C(j-1,:)=0;
n=n-1;
case 5 % Customer defined boundary condition. Must be linear.
if isnumeric(bd{i,4})
m=bd{i,4};
else
m=bd{i,4}(h);
end
p=length(m);
A(l-1+n,j-p+1:j)=m;
case 6
if isnumeric(bd{i,4})
m=bd{i,4};
else
m=bd{i,4}(h);
end
p=length(m);
A(l-1+n,j:j+p-1)=m;
case 7 %
if isnumeric(bd{i,4})
m=bd{i,4};
else
m=bd{i,4}(h);
end % length of m must be odd
p=(length(m)-1)/2;
A(l-1+n,j-p:j+p)=m;
end
n=n+1;
end
if n>2
fprintf('Too many boundary conditions!')
end
% For test
%A=full(A);B=full(B);C=full(C);
end
u2=A\(B*u1+C);
while n_tout<=length(tlist) && (tlist(n_tout)>=t(tstep) && tlist(n_tout)<=t(tstep+1))
% Algorithm to find values for each x on xlist
n_xout=1; % outputing the N0. n_xout value.
for i=1:l-1 % Search in all x intervals.
while n_xout<=length(xlist) && (x(i)<=xlist(n_xout) && x(i+1)>=xlist(n_xout))
% desired x value between x(i) and x(i+1)
uout1(1,n_xout)=u1(i+1:-1:i)'*[xlist(n_xout)-x(i);x(i+1)-xlist(n_xout)]/(x(i+1)-x(i));
n_xout=n_xout+1;
end
if n_xout>length(xlist)
break;
end
end
% Algorithm to find values for each x on xlist
n_xout=1; % outputing the N0. n_xout value.
for i=1:l-1 % Search in all x intervals.
while n_xout<=length(xlist) && (x(i)<=xlist(n_xout) && x(i+1)>=xlist(n_xout))
% desired x value between x(i) and x(i+1)
uout2(1,n_xout)=u2(i+1:-1:i)'*[xlist(n_xout)-x(i);x(i+1)-xlist(n_xout)]/(x(i+1)-x(i));
n_xout=n_xout+1;
end
if n_xout>length(xlist)
break;
end
end
uout(n_tout,:)=[tlist(n_tout)-t(tstep),t(tstep+1)-tlist(n_tout)]/(t(tstep+1)-t(tstep))*[uout2;uout1];
n_tout=n_tout+1;
end
if n_tout>length(tlist)
break
end
u1=u2;
end

