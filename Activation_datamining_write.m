%% This file is only appropriate for TGA_redox cycles activation data mining
% This file is only for data saved as xlsx files
% Only be available when the oxidation rate of sample is good enough
% Only be available when iron oxide as active materials and oxygen carrier
% Change the conversion calculation when used as other materials
% The 100% conversion means Fe3O4
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
t_red=15;% reduction duration (min)
t_flush=5; % flush duration (min)
t_ox=5; % oxidation duration (min)

%data read from excel
[FileName,PathName,~]=uigetfile('*.xls*','open an xlsx document saved as xlsx');
data=xlsread([PathName,FileName]);
t=data(:,2)*60; % Based on min
temp=data(:,3); % sample temperature
dtg=data(:,5);
tg=data(:,4); % can be used to calculate conversion
%sample name create
 Samplename=input('What is the name of your sample?...','s');
 temperature=input('What is the temperature of your reaction_C?...');
%data for conversion calculation
 dopant=input('What is the dopant: non=0, La=1, Co=2, Ni=3, Cu=4?');
 percent=input('What is the percentage of dopant?');

%Judge_1 if the number is okay for thres
judge_1=0;
%Judge_2 if the number for standard and dtg is good
judge_2(1:n_cycle)=0;
while judge_1<0.5;
    
% findpeaks(dtg,'MinPeakProminence',thres) % for testing
thres=0.01*input('What is the minimum peak prominence in dTG? 1 means 1e-2?...');
% minimum peak prominence in dTG
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

I_ox2=I_p+dI_ox;
tg_ox=tg(I_ox2);
I_re=I_p-dI_re;
tg_re=tg(I_re);
I_oxbegin=I_p(1)-dI_p+dI_ox;
tg_oxbegin=[tg(I_oxbegin);tg_ox(1:end-1)];

%Figure plot to determine the places is write or not for conversion
%calculation

figure
plot (t,dtg,...
    t(I_p),dtg(I_p),'*k')

figure
plot(t,tg,...
    t([I_oxbegin;I_ox2]),tg([I_oxbegin;I_ox2]),'*r',...
    t(I_re),tg(I_re),'*b')

%Judge if it is okay or not and then run next step
judge_1=input('Input 1 if it is okay for points, Input 0 if it is not..........');
end

I_ox=[I_oxbegin;I_ox2(1:(n_cycle-1))];
for (i=1:n_cycle)

while judge_2(i)<0.5;
%depends on min........
%Find all this cycle value
tg_1cycle=tg(I_ox(i):I_re(i));

%Based on min
Time=t(I_ox(i):I_re(i));
dtg_1cycle=dtg(I_ox(i):I_re(i));

%This startindex value is to find appropriate starting point of TGA data
startindex(i)=input(['the start indicate for dtg?...'...
    num2str(i),'_cycles...']);

%Set the beginning for the TGA curve
[~,I_max(i)]=max(tg_1cycle);

%Use the creteria to find the point again
[I_redstart1(i),~]=find(dtg_1cycle(I_max(i):end)<startindex(i),1);

%Obtain all the point
I_redstart(i)=I_redstart1(i)+I_max(i)-1;
%Get the meaningful I_redstop

I_redstop(i)=I_redstart(i)+round(t_red/dt)-1;%based on min

%Make sure the I_redstop doesn't exceed
if (I_redstop(i)>length(Time))
    I_redstop(i)=length(Time);
end
%Plot the data acquisition plot to see if it is okay
figure
subplot(2,1,1)
plot(Time,tg_1cycle,...
    Time(I_redstart(i):I_redstop(i)),tg_1cycle(I_redstart(i):I_redstop(i)),'*r');
subplot(2,1,2)
plot(Time,dtg_1cycle,...
    Time(I_redstart(i)),dtg_1cycle(I_redstart(i)),'*r');

standard(i)=input(['Input the standard of your choice to calibrate the m and k...',...
    num2str(i),'_cycles...']);
%Calculate the m and k for this line, conversion is based on Fe3O4
TGAdata=[Time(I_redstart(i):I_redstop(i)),tg_1cycle(I_redstart(i):I_redstop(i))];
[m(i,:),k(i,:),mr2(i,:),kR2(i,:)] = kcalculate(TGAdata,dopant,percent,standard(i));
judge_2(i)=input('Input 1 if it is okay for points, Input 0 if it is not..........');
end

% Write all the calculation results back to the excel files
% construct the cell for conversion information
Conclusion=num2cell([(1:i)',m,k,mr2,kR2]);
% Results=num2cell(Conclusion);
Title_Results={'Cycle_number','m','k_m=0.54','k_m=0.57',...
    'k_m=0.62','k_m=1','k_m=1.07','k_m=1.11','k_m=2',...
    'k_m=3','mr2','kR2_m=0.54','kR2_m=0.57',...
    'kR2_m=0.62','kR2_m=1','kR2_m=1.07','kR2_m=1.11','kR2_m=2',...
    'kR2_m=3'};
Datatoexcel=[Title_Results;Conclusion];

%Write all the data to the first conversion file

xlswrite([PathName,Samplename],Datatoexcel,num2str(temperature),'A1');

end

