function Cm = MCco()
%MCCO Summary of this function goes here
%generate the monochromatic coefficient matrix
%the coefficient take into account the scattering factor and the
%integration weights
%   Detailed explanation goes here
m=12;
IA=linspace(1,0,m/2+2);
IA=IA(:,2:(m/2+1));
um=[IA,flip(-IA)];
x1=IA;
y1=11.53827*IA+18.91319;
x2=-IA;
y2=-11.53827*IA+18.91319;
x3=-0.118*log(IA)-0.2568;
y3=-15.836*IA+16.238;
Cms=zeros(m,m);
%row represents a specific scattering intensity from several incident angle
%column represents a specific incident angle to several scattering angle
for i=1:m/2
    if i==1
        X=[x1(i)^2,x1(i),1;x2(i)^2,x2(i),1;x3(i)^2,x3(i),1];
        A=X\[y1(i);y2(i);y3(i)];
        Cms(:,i)=transpose(A)*[um.^2;um;ones(1,m)];
        Cms(:,m+1-i)=flip(Cms(:,i));
    else if i==m/2
            Ap=y1(i)/log(x1(i));
            An=y2(i)/log(-x2(i));
            Cms(:,i)=[Ap*log(IA),An*log(flip(IA))];
            Cms(:,m+1-i)=flip(Cms(:,i));
        else
            X=[x1(i)^2,x1(i),1;x2(i)^2,x2(i),1;x3(i)^2,x3(i),1];
            A=X\[y1(i);y2(i);y3(i)];
            Ap=y1(i)/log(x1(i));
            An=y2(i)/log(-x2(i));
            Cms(:,i)=[Ap*log(IA(1:i)),transpose(A)*[um(:,(i+1):(m-i)).^2;um(:,(i+1):(m-i));ones(1,m-2*i)],An*log(flip(IA(1:i)))];
            Cms(:,m+1-i)=flip(Cms(:,i));
        end
    end
end

QWs=NCQWO5();
Cm=Cms*diag([QWs,QWs]);

end