%% This file is only appropriate for TGA_redox cycles conversion calculation

clear
clc
close('all')

%data read from excel
[FileName,PathName,~]=uigetfile('*.xls*','select the file directly from TGA');
data=readtable([PathName,FileName],'FileType','text','Delimiter','\t','HeaderLines',12);
t=table2array(data(:,2)); % time (h)
temp=table2array(data(:,3)); % sample temperature
dtg=table2array(data(:,5));
tg=table2array(data(:,4)); % can be used to calculate conversion

%data for conversion calculation
Support=input('what is the weight percentage(wt%) of support?');
initial_weight=input('what is the initial weight of your sample in TGA(mg)?');
high_T_weight=input('what is the weight of your sample at higher temperature in TGA(mg)?')

[weightper,c,index]=TGA_convcalc_gen1particle(tg,Support,initial_weight,high_T_weight); % conversion calculation
%Set the time of beginning
t=t-t(index);

figure
plot(t,c)

%Write all the calculation results back to the excel files
newFileName=strrep(FileName,FileName((end-3):end),'_forkinetics.xlsx');

%construct the cell for conversion information
Conversion=table(t,tg,dtg,c);
Conversion_cell=table2cell(Conversion);
Title_Conversion={'Time_(min)','TG_(mg)',...
    'DTG_(mg/min)','Conversion_(%)'};
Datatoexcel_1=[Title_Conversion;Conversion_cell];

%Write all the data to the first conversion file
xlswrite([PathName,newFileName],Datatoexcel_1,'Conversion','A1');

%Clean conversion file
T=table(t(index:end),c(index:end));
cleanconversion=table2cell(T);
Title_cleanconversion={'Time(min)','Conversion(%)'};
Datatoexcel_2=[Title_cleanconversion;cleanconversion];
%Write redox cycles to another file
xlswrite([PathName,newFileName],Datatoexcel_2,'Clean_conversion','A1');