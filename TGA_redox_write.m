%% This file is only appropriate for TGA_redox cycles conversion calculation
% This file is only used for xlsx document, which is not directly from TGA software 
% Only be available when the oxidation rate of sample is good enough
% Only be available when iron oxide as active materials and oxygen carrier
% Change the conversion calculation when used as other materials
% Before using this program, check the following status of your TGA file
% First, if you have the DTG data; Second, if the time is in the hour unit?
% If not, first: get the DTG data from the TGA instrument;
% If not, second: delete the 60 for min or change it to *60 for seconds
% dI_ox=round((t_ox+t_flush/2)/60/dt);
% dI_re=round((t_flush/2)/60/dt); 
% afterward change the output line
% Third, check cycles,thres,flush time and oxidation time of each
% experiments and revise the parameters section below

clear
clc
close('all')


% Parameters
n_cycle=15; % number of redox cycles
thres=1e-2; % minimum peak prominence in dTG
t_flush=5; % flush duration (min)
t_ox=5; % oxidation duration (min)

%data read from excel
[FileName,PathName,~]=uigetfile('*.xls*','open an xlsx document not directly from TGA software');
data=xlsread([PathName,FileName]);
t=data(:,2); % time (h)
temp=data(:,3); % sample temperature
dtg=data(:,5);
tg=data(:,4); % can be used to calculate conversion

%data for conversion calculation
dopant=input('what is the dopant: non=0, La=1, Co=2, Ni=3, Cu=4?');
percent=input('what is the percentage of dopant?');
[weightper,c]=TGA_convcalc(tg,dopant,percent); % conversion calculation

% findpeaks(dtg,'MinPeakProminence',thres) % for testing
[pks,I_p]=findpeaks(dtg,'MinPeakProminence',thres);

% find the reducing and oxidation conversion of each cycles
% using find peaks and index to determine the places
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

%Figure plot to determine the places is write or not for conversion
%calculation
figure
plot(t,[c,tg,dtg],...
    t(I_p),dtg(I_p),'*k',...
    t([I_ox2;I_ox]),c([I_ox2;I_ox]),'*r',...
    t(I_re),c(I_re),'*b')

figure
bar([RE,OX])
xlabel('Cycle'); ylabel('Conversion (%)');
legend('Reduction','Oxidation')

%Write all the calculation results back to the excel files
Conversion=table(t,temp,tg,dtg,weightper,c);
Conversion_cell=table2cell(Conversion);
Title_Conversion={'Time_(h)','Temperature_(C)','TG_(mg)',...
    'DTG_(mg/min)','Weightper_(%)','Conversion_(%)'};
Datatoexcel_1=[Title_Conversion;Conversion_cell];

%Write all the data to the first conversion file
xlswrite([PathName,FileName],Datatoexcel_1,'Conversion','A1');

%Redox cycles reduction and oxidation
T=table([1:n_cycle]',RE,OX);
Redox_cycles_cell=table2cell(T);
Title_Redox_cycles={'Cycles','Reduction(%)','Oxidation(%)'};
Datatoexcel_2=[Title_Redox_cycles;Redox_cycles_cell];
%Write redox cycles to another file
xlswrite([PathName,FileName],Datatoexcel_2,'Redox_cycles','A1');