%% This file is only appropriate for TGA_redox cycles conversion calculation
% This file is only for data directly from TGA software 
% Only be available when the oxidation rate of sample is good enough
% Only be available when iron oxide and copper oxide as active materials and oxygen carrier
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
n_cycle=10; % number of redox cycles
thres=1e-2; % minimum peak prominence in dTG
t_flush=5; % flush duration (min)
t_ox=15; % oxidation duration (min)

t_flush=t_flush/60;
t_ox=t_ox/60;

%data read from excel
[FileName,PathName,~]=uigetfile('*.xls*','select the file directly from TGA');
data=readtable([PathName,FileName],'FileType','text','Delimiter','\t','HeaderLines',12);
t=table2array(data(:,2)); % time (h)
temp=table2array(data(:,3)); % sample temperature
dtg=table2array(data(:,5));
tg=table2array(data(:,4)); % can be used to calculate conversion

%data for conversion calculation
Fe2O3=input('what is the weight percentage of Fe2O3(%)?...');
CuO=input('what is the weight percentage of CuO(%)?...');

[weightper,c] = TGA_convcalc_CuO_Al2O3_Fe2O3(tg,Fe2O3,CuO); % conversion calculation

% findpeaks(dtg,'MinPeakProminence',thres) % for testing
[pks,I_p]=findpeaks(dtg,'MinPeakProminence',thres);

% find the reducing and oxidation conversion of each cycles
% using find peaks and index to determine the places
[~,temp_I]=sort(pks,'descend');
temp_I2=temp_I(n_cycle+1:end); % buggy when length(temp_I) = n_cycle
I_p(temp_I2)=[];

dt=mean(diff(t));
dI_p=mode(diff(I_p));
dI_ox=round((t_ox+t_flush/2)/dt);
dI_re=round((t_flush/2)/dt);
I_pmax=I_p-round(t_flush/dt); % determine the beginning of flush

% Linear regression
I_oxmax=I_p-round((t_flush/2+1/60)/dt); % one minute ahead of the middle time
I_oxmin=I_p-round((t_flush/2-1/60)/dt); % one minute after the middle time, use the two minute data to do linear regression

for i=1:n_cycle
    I=I_oxmax(i):I_oxmin(i);
    I=I';
    c_ox=c(I);
    X=[ones(size(I)),I];
    [b,bint,r,rint,stat]=regress(c_ox,X); % Linear regression
    Iini=I_pmax(i):I_p(i);
    Iini=Iini';
    cini=c(Iini);
    [c_oxini(i,:),I_oxini(i,:)]=max(cini); % Find the initial flushing index
    I_oxini(i,:)=I_oxini(i,:)+I_pmax(i)-1; % Find the initial flushing index
    c_re(i,:)=b(1)+b(2)*I_oxini(i); % Find the conversion at the initial flushing time
    R(i)=stat(1)^0.5; % Variance for linear regression
end

I_ox=I_p+dI_ox;
I_re=I_oxini;
c_ox=c(I_ox);
I_ox2=I_p(1)-dI_p+dI_ox;
c_ox2=[c(I_ox2);c_ox(1:end-1)];
OX=c_re-c_ox;
RE=c_re-c_ox2;
R=R';
%Figure plot to determine the places are right or not for conversion
%calculation
figure
plot(t,[c,tg,dtg],...
    t(I_p),dtg(I_p),'*k',...
    t([I_ox2;I_ox]),c([I_ox2;I_ox]),'*r',...
    t(I_re),c_re,'*b')

figure
bar([RE,OX])
xlabel('Cycle'); ylabel('Conversion (%)');
legend('Reduction','Oxidation')

%Write all the calculation results back to the excel files
newFileName=strrep(FileName,FileName((end-3):end),'_handled_with oxygen in dopants.xlsx');

%construct the cell for conversion information
Conversion=table(t,temp,tg,dtg,weightper,c);
Conversion_cell=table2cell(Conversion);
Title_Conversion={'Time_(h)','Temperature_(C)','TG_(mg)',...
    'DTG_(mg/min)','Weightper_(%)','Conversion_(%)'};
Datatoexcel_1=[Title_Conversion;Conversion_cell];

%Write all the data to the first conversion file
xlswrite([PathName,newFileName],Datatoexcel_1,'Conversion','A1');

%Redox cycles reduction and oxidation
T=table([1:n_cycle]',RE,OX);
Redox_cycles_cell=table2cell(T);
Title_Redox_cycles={'Cycles','Reduction(%)','Oxidation(%)'};
Datatoexcel_2=[Title_Redox_cycles;Redox_cycles_cell];
%Write redox cycles to another file
xlswrite([PathName,newFileName],Datatoexcel_2,'Redox_cycles','A1');