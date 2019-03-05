% Solve du/dt = d^2u/dx^2 + d^2u/dy^2 on a rectangular domain with ADI
% (Alternating Direction Implicit) method
X = [0,1];
Y = [0,1];
T = [0,0.01];
u_sln = @(x,y,t)sin(pi*x).*sin(3*pi*y)*exp(-10*pi^2*t);
u0 = @(x,y)sin(pi*x).*sin(3*pi*y); %initial conditon
ub = @(x,y,t)0;
klist = (T(2)-T(1))*0.5.^(0:9);
mlist = (10*2.^(0:7));
hlist = (X(2)-X(1))./mlist;
err = zeros(length(klist),length(mlist));
for i = 1:length(klist)
k = klist(i);
for j = 1:length(mlist)
m = mlist(j);
n = m;
hx = (X(2)-X(1))/(m+1);
x = X(1):hx:X(2);
hy = (Y(2)-Y(1))/(n+1);
y = Y(1):hy:Y(2);
[xgrid,ygrid] = meshgrid(x,y);
Ax = gallery('tridiag',m,k/2/hx^2,-k/hx^2,k/2/hx^2);
Ay = gallery('tridiag',n,k/2/hy^2,-k/hy^2,k/2/hy^2);
Im = speye(m);
In = speye(n);
% This step will make the first index represent x, and second index represent y.
xgrid = xgrid';
ygrid = ygrid';
u1 = u0(xgrid,ygrid);
v1 = u1(2:m+1,2:n+1);
t = T(1);
n_step = 2;
% uout = zeros((m+2)*(n+2),ceil((T(2)-T(1))/k)+1);
% uout(:,1) = u1(:);
while t < T(2)
Bx = k/2/hx^2*(sparse(1,1:n,ub(xgrid(1,:),ygrid(1,:),t+k/2),m,n)+...
sparse(m,1:n,ub(xgrid(m,:),ygrid(m,:),t+k/2),m,n));
By1 = k/2/hy^2*(sparse(1:m,1,ub(xgrid(:,1),ygrid(:,1),t),m,n)+...
sparse(1:m,n,ub(xgrid(:,n),ygrid(:,n),t),m,n));
By2 = k/2/hy^2*(sparse(1:m,1,ub(xgrid(:,1),ygrid(:,1),t+k),m,n)+...
sparse(1:m,n,ub(xgrid(:,n),ygrid(:,n),t+k),m,n));
v2 = ((Im+Ax)*((Im-Ax)\(v1*(In+Ay)+By1+Bx))+Bx+By2)/(In-Ay);
v1 = v2;
u1(2:m+1,2:n+1) = v2;
u1(1,:) = ub(xgrid(1,:),ygrid(1,:),t+k);
u1(m+2,:) = ub(xgrid(m,:),ygrid(m,:),t+k);
u1(:,1) = ub(xgrid(:,1),ygrid(:,1),t+k);
u1(:,n+2) = ub(xgrid(:,n),ygrid(:,n),t+k);
% uout(:,n_step) = u1(:);
t = t+k;
n_step = n_step+1;
end
% for i = 1:size(uout,2)
% mesh(xgrid,ygrid,reshape(uout(:,i),m+2,n+2))
% pause(0.1)
% end
err1 = u1-u_sln(xgrid,ygrid,t);
err(i,j) = max(abs(err1(:)));
end
end
[K,H]=meshgrid(klist,hlist);
subplot(1,2,1)
surf(K,H,err');
set(gca,'XScale','log','YScale','log','ZScale','log')
xlabel('k')
ylabel('h')
zlabel('error')
err1 = zeros(1,length(hlist));
for i=1:min(length(hlist),length(klist))
err1(i)=err(i,i);
end
p = polyfit(log(hlist(4:8)),log(abs(err1(4:8))),1);
subplot(1,2,2)
loglog(hlist,err1,'o-')
ylabel('Error')
xlabel('h')