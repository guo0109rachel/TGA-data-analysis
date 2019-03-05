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

%smooth H2 concentration
MSgasraw(:,5)=smooth(MSgasraw(:,5),30,'moving');

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
TGA_MSresults2(:,10:12)=TGA_MSresults2(:,10:12)/WeightO*1000;
%Change mg/min to %(O)/min
TGA_MSresults2(:,6)=TGA_MSresults2(:,6)/WeightO*100;

%Change initial concentration to 0
initialC=TGA_MSresults2(1:2100,10:12);%find the average concentration of first 400C
averageC=mean(initialC,1);
averageC=[zeros(1,9),averageC];
TGA_MSresults2=TGA_MSresults2-averageC;%substract the average concentration of first 400C
TGA_MSresults2=array2table(TGA_MSresults2);

%write it to the file
 HEAD=[{'time_s','time_min','TGA_temp_C','TGA_weight_mg','TGA_weightper',...
   'TGA_DTG','Conversion'},MSfinalhead];
TGA_MSdata=table2cell(TGA_MSresults);
Finalresults=[HEAD;TGA_MSdata];
TGA_MSdata2=table2cell(TGA_MSresults2);
Finalresults2=[HEAD;TGA_MSdata2];

%write it to a new excel file
newFileName=strrep(TGAFileName,TGAFileName((end-3):end),'_smoothH30.xlsx');
xlswrite([TGAPathName,newFileName],Finalresults,'original data_handled','A1');
xlswrite([TGAPathName,newFileName],Finalresults2,'modified data_handled','A1');