%% 
clear
clc
close('all')

% Parameters
n_cycle=15; % number of redox cycles
thres=1e-2; % minimum peak prominence in dTG
t_flush=5; % flush duration (min)
t_ox=5; % oxidation duration (min)

[FileName,PathName,~]=uigetfile('*.xls*');
data=xlsread([PathName,FileName]);
t=data(:,2); % time (h)
c=data(:,7); % conversion
dtg=data(:,5);
tg=data(:,4); % can be used to calculate conversion

% findpeaks(dtg,'MinPeakProminence',thres) % for testing

[pks,I_p]=findpeaks(dtg,'MinPeakProminence',thres);

[~,temp_I]=sort(pks,'descend');
temp_I2=temp_I(n_cycle+1:end); % buggy when length(temp_I) = n_cycle
I_p(temp_I2)=[];

dt=mean(diff(t));
dI_p=mode(diff(I_p));
dI_ox=round((t_ox+t_flush/2)/60/dt);
dI_re=round((t_flush/2)/60/dt);

I_ox=I_p+dI_ox;
c_ox=c(I_ox);
I_re=I_p-dI_re;
c_re=c(I_re);
I_ox2=I_p(1)-dI_p+dI_ox;
c_ox2=[c(I_ox2);c_ox(1:end-1)];
OX=c_re-c_ox;
RE=c_re-c_ox2;



figure
plot(t,[c,tg,dtg],...
    t(I_p),dtg(I_p),'*k',...
    t([I_ox2;I_ox]),c([I_ox2;I_ox]),'*r',...
    t(I_re),c(I_re),'*b')

figure
bar([RE,OX])
xlabel('Cycle'); ylabel('Conversion (%)');
legend('Reduction','Oxidation')


