%% This file is only appropriate for TGA_MS TPR Operation
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
[TGAFileName,TGAPathName,~]=uigetfile('*.xls*','select the TGA data for the same one with MS');
TGAdata=xlsread([TGAPathName,TGAFileName]);

%get the useful data from MS
%get the timeline of MS
tms=strrep(tms,'"','');
tms_raw=datenum(tms,'mm/dd/yyyy HH:MM:SS PM');
TMS=(tms_raw-tms_raw(1)).*3600.*24;%time unit is second

%get the conversion line
MShead=strrep(MShead,'"','');
tf=contains(MShead,'%');
[~,column]=find(tf);
%CO2 is the first one, Helium is the second, CO and H2 and then CH4; change
%data in TGA_MSfunc_direct program;
MSgasraw=MSdata(:,(column(1:end)-1));
%Construct the head of the MS part
MSfinalhead=MShead(column(1:end));

%TGA and MS data combination
dopant=input('what is the dopant: non=0, La=1, Co=2, Ni=3, Cu=4?......');
percent=input('what is the percentage of dopant?......');
delay=input('How long did you wait for turing on MS after TGA(s)?......');
Helium_flow=input('What is the flow rate of Helium? (mL/min)......');
TGA_MSresults = TGA_MSfunc_direct(TGAdata,TMS,MSgasraw,dopant,percent,Helium_flow,delay);

%write it to the file
 HEAD=[{'time_s','time_min','TGA_temp_C','TGA_weight_mg','TGA_weightper',...
   'TGA_DTG','Conversion'},MSfinalhead,{'CO2_flow','CO_flow','H2_flow','CH4_flow',...
   'CO2_percent','CO_percent','H2_percent','CH4_percent'}];
TGA_MSdata=table2cell(TGA_MSresults);
Finalresults=[HEAD;TGA_MSdata];

%write it to the excel file
xlswrite([TGAPathName,TGAFileName],Finalresults,'data_handled','A1');



