function [m,k,mr2,kR2] = kcalculate(TGAdata,dopant,percent,standard)
%Calculate the k of the reaction for specific temperature

%Get the time from from TGAdata
%Notice that time is depend on minute basis
Timeline=(TGAdata(:,1)+TGAdata(2,1)-TGAdata(1,1)-TGAdata(1,1));

%Calculate the conversion based on reduced to Fe3O4
rawdata=TGAdata(:,2);
Conversion=TGA_conv_Fe3O4(rawdata,dopant,percent);
Conversion_index=[Conversion;1];

%Delete the conversion data greater or equal to standard
i=find(Conversion_index>=standard,1);
Conversion(i:end)=[];
t1=Timeline(2:(i-1));
X1=Conversion(2:end);
n1=length(t1);
%polyfit to calculate m and r2
[c,~,~,~,s]=regress(log(-log(1-X1)),[ones(n1,1),log(t1)]);
m=c(2);
mr2=s(1);

%Fitting for k
X2=Conversion;
t2=Timeline(1:(i-1));
n2=length(t2);

%when m=0.54 get k(1) for three-dimensional diffusion
[b1,~,~,~,s1]=regress((1-(1-X2).^(1/3)).^2,[ones(n2,1),t2]);
k(1)=b1(2);
kR2(1)=s1(1);

%when m=0.57 get k(2) for two-dimensional diffusion
[b2,~,~,~,s2]=regress(-(1-X2).*log(1-X2),[ones(n2,1),t2]);
k(2)=b2(2);
kR2(2)=s2(1);

%when m=0.62 get k(3) for one-dimensional diffusion
[b3,~,~,~,s3]=regress(X2.^2,[ones(n2,1),t2]);
k(3)=b3(2);
kR2(3)=s3(1);

%when m=1 get k(4) for first order reaction;
[b4,~,~,~,s4]=regress(-log(1-X2),[ones(n2,1),t2]);
k(4)=b4(2);
kR2(4)=s4(1);

%when m=1.07 get k(5) for sphere Phase-boundary Controlled
[b5,~,~,~,s5]=regress(1-(1-X2).^(1/3),[ones(n2,1),t2]);
k(5)=b5(2);
kR2(5)=s5(1);

%when m=1.11 get k(6) for cylinder Phase-boundary Controlled
[b6,~,~,~,s6]=regress(1-(1-X2).^0.5,[ones(n2,1),t2]);
k(6)=b6(2);
kR2(6)=s6(1);

%when m=2 get k(7) for two-dimensional growth of nuclei
[b7,~,~,~,s7]=regress((-log(1-X2)).^0.5,[ones(n2,1),t2]);
k(7)=b7(2);
kR2(7)=s7(1);

%when m=3 get k(8) for three-dimensional growth of nuclei
[b8,~,~,~,s8]=regress((-log(1-X2)).^(1/3),[ones(n2,1),t2]);
k(8)=b8(2);
kR2(8)=s8(1);


%plotfigure
figure
%Plot for m results
subplot(3,3,1);
plot(log(t1),m.*log(t1)+c(1),...
    log(t1),log(-log(1-X1)),'*g');

title(['fitting m=',num2str(m),',r2=',num2str(mr2)]);

%Plot for k results when m=0.54 get for three-dimensional diffusion 
subplot(3,3,2);
plot(t2,(1-(1-X2).^(1/3)).^2,'*g',...
t2,k(1).*t2+b1(1));
title(['m=0.54,k=',num2str(k(1)),...
    ',R2=',num2str(kR2(1))]);

%Plot for k results when m=0.57 get for two-dimensional diffusion
subplot(3,3,3);
plot(t2,-(1-X2).*log(1-X2),'*g',...
    t2,k(2).*t2+b2(1));
title(['m=0.57,k=',num2str(k(2)),...
    ',R2=',num2str(kR2(2))]);

%Plot for k results when m=0.62 get for one-dimensional diffusion
subplot(3,3,4);
plot(t2,X2.^2,'*g',...
    t2,k(3).*t2+b3(1));
title(['m=0.62,k=',num2str(k(3)),...
    ',R2=',num2str(kR2(3))]);

%Plot for k results when m=1 for first order reaction;
subplot(3,3,5);
plot(t2,-log(1-X2),'*g',...
    t2,k(4).*t2+b4(1));
title(['m=1,k=',num2str(k(4)),...
    ',R2=',num2str(kR2(4))]);

%Plot for k results when m=1.07 for sphere Phase-boundary Controlled
subplot(3,3,6);
plot(t2,1-(1-X2).^(1/3),'*g',...
    t2,k(5).*t2+b5(1));
title(['m=1.07,k=',num2str(k(5)),...
    ',R2=',num2str(kR2(5))]);

%Plot for k results when m=1.11 for cylinder Phase-boundary Controlled
subplot(3,3,7);
plot(t2,1-(1-X2).^0.5,'*g',...
    t2,k(6).*t2+b6(1));
title(['m=1.11,k=',num2str(k(6)),...
    ',R2=',num2str(kR2(6))]);

%Plot for k results when m=2 for two-dimensional growth of nuclei
subplot(3,3,8);
plot(t2,(-log(1-X2)).^0.5,'*g',...
    t2,k(7).*t2+b7(1));
title(['m=2,k=',num2str(k(7)),...
    ',R2=',num2str(kR2(7))]);

%Plot for k results when m=3 for three-dimensional growth of nuclei
subplot(3,3,9);
plot(t2,(-log(1-X2)).^(1/3),'*g',...
    t2,k(8).*t2+b8(1));
title(['m=3,k=',num2str(k(8)),...
    ',R2=',num2str(kR2(8))]);
end