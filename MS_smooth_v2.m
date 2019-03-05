%% This file is only appropriate for TGA_MS TPR Operation
% Read the file directly from TGA software
%Instrument
clear
clc
close('all')
%data read from MS-txt
[MSFileName,MSPathName,~]=uigetfile('*.txt*','select the MS data');
MSfile=importdata([MSPathName,MSFileName],' ');
MSfile(1:53)=[]; MSfile(end)=[];
MShead=strsplit(MSfile{1},'\t');

%Create a datafile to read only numbers
MSdat=cell(length(MSfile)-1,length(MShead));
for i=1:length(MSfile)-1
    MSdat(i,:)=strsplit(MSfile{i+1},'\t');
end
tms=MSdat(:,1);
MSdata=str2double(MSdat(:,2:end));

%data read from TGA
[TGAFileName,TGAPathName,~]=uigetfile('*.xls*','select the file directly from TGA software');
readingdata=readtable([TGAPathName,TGAFileName],'FileType','text','Delimiter','\t','HeaderLines',11);
TGAdata=table2array(readingdata);

%get the useful data from MS
%get the timeline of MS
tms=strrep(tms,'"','');
tms_raw=datenum(tms,'mm/dd/yyyy HH:MM:SS PM');
TMS=(tms_raw-tms_raw(1)).*3600.*24;%time unit is second

%get the conversion line
MShead=strrep(MShead,'"','');
tf=contains(MShead,'%');
[~,column]=find(tf);
MSgasraw=MSdata(:,(column(1:end)-1));
%Construct the head of the MS part
MSfinalhead=[MShead(column(1:end))];

%smooth concentration

%TGA and MS data combination
dopant=input('what is the dopant: non=0, La=1, Co=2, Ni=3, Cu=4?');
percent=input('what is the percentage of dopant?');
delay=input('How long did you wait for turing on MS after TGA(s)?');
TGA_MSresults = TGA_MSfunc_direct(TGAdata,TMS,MSgasraw,dopant,percent,delay);

%Calculate total mass of O
tg=TGAdata(:,4);
WeightO=WeightO(tg,dopant,percent);

%Change concentration to %(gO)/min
TGA_MSresults2=table2array(TGA_MSresults);
temp = table2array(TGA_MSresults(:,3));
TGA_MSresults2(:,10:12)=TGA_MSresults2(:,10:12)/WeightO*1000;
H2 = TGA_MSresults2(:,10);
CO = TGA_MSresults2(:,11);
CO2 = TGA_MSresults2(:,12);

n = 0;
while n ~=1
    k = input('parameter to smooth H2     k =  ');
    H2f  = smooth(temp, H2, k, 'lowess');
    plot(temp, H2), hold on
    plot(temp, H2f, 'r'), hold off
    n = input('if the result is satisfactory, enter 1      ');
end

n = 0;
while n ~=1
    k = input('parameter to smooth CO     k =  ');
    COf  = smooth(temp, CO, k, 'lowess');
    plot(temp, CO), hold on
    plot(temp,COf, 'r'), hold off
    n = input('if the result is satisfactory, enter 1      ');
end

n = 0;
while n ~=1
    k = input('parameter to smooth CO2     k =  ');
    CO2f  = smooth(temp, CO2, k, 'lowess');
    plot(temp, CO2), hold on
    plot(temp, CO2f, 'r'), hold off
    n = input('if the result is satisfactory, enter 1      ');
end


n = 0;
while n ~=1
    a = input('the numer of the original point you want to plot   ')
    subplot(3,1,1)
    plot(temp,H2f), hold on
    plot(temp(1:a:end),H2(1:a:end),'*'), hold off
    title('H2 concentration')
    subplot(3,1,2)
    plot(temp,COf), hold on
    plot(temp(1:a:end),CO(1:a:end),'*'), hold off
    title('CO concentration')
    subplot(3,1,3)
    plot(temp,CO2f), hold on
    plot(temp(1:a:end),CO2(1:a:end),'*'), hold off
    title('CO2 concentration')
    n = input('if the result is satisfactory, enter 1      ');
end

% change the initial concentration to 0
initialC=TGA_MSresults2(1:1050,10:12);%find the average concentration of first 400C
averageC1=mean(initialC,1);
averageC=[zeros(1,9),averageC1];
TGA_MSresults2=TGA_MSresults2-averageC;%substract the average concentration of first 400C

% save the original data to plot
H2 = TGA_MSresults2(:,10);
CO = TGA_MSresults2(:,11);
CO2 = TGA_MSresults2(:,12);

%Change mg/min to %(O)/min
TGA_MSresults2(:,6)=TGA_MSresults2(:,6)/WeightO*100;
TGA_MSresults4=TGA_MSresults2;%substract the average concentration of first 400C
TGA_MSresults4(:,10:12) = [H2f, COf, CO2f];
C2=min(TGA_MSresults4(:,10:12));
TGA_MSresults4(:,10:12) = [H2f, COf, CO2f] - C2;% change initial concentration to zero
TGA_MSresults4 = array2table(TGA_MSresults4);
TGA_MSresults2=array2table(TGA_MSresults2);

% data to plot scatter
l_temp = length(temp);
i_temp = 1:a:l_temp;
TGA_MSresults3 = [temp(i_temp),H2(i_temp),CO(i_temp),CO2(i_temp)];
TGA_MSresults3 = array2table(TGA_MSresults3);

%write it to the file
 HEAD=[{'time_s','time_min','TGA_temp_C','TGA_weight_mg','TGA_weightper',...
   'TGA_DTG','Conversion'},MSfinalhead];
HEAD3 = [{'TGA_temp_C'},MSfinalhead(3:5)];
TGA_MSdata=table2cell(TGA_MSresults);
Finalresults=[HEAD;TGA_MSdata];
TGA_MSdata2=table2cell(TGA_MSresults2);
Finalresults2=[HEAD;TGA_MSdata2];
TGA_MSdata3 = table2cell(TGA_MSresults3);
Finalresults3 = [HEAD3;TGA_MSdata3];
TGA_MSdata4 = table2cell(TGA_MSresults4);
Finalresults4 = [HEAD;TGA_MSdata4];

%write it to a new excel file
newFileName=strrep(TGAFileName,TGAFileName((end-3):end),'_smooth v2.xlsx');
xlswrite([TGAPathName,newFileName],Finalresults,'original data_handled','A1');
xlswrite([TGAPathName,newFileName],Finalresults2,'modified data_handled','A1');
xlswrite([TGAPathName,newFileName],Finalresults4,'smooth v2_handled','A1');
xlswrite([TGAPathName,newFileName],Finalresults3,'points to plot','A1');