%%
clear
clc
close('all')

[FileName,PathName,~]=uigetfile('*.xls');
dat=readtable([PathName,FileName],'FileType','text','Delimiter','\t','HeaderLines',11);
t=dat{:,2};
tg=dat{:,4};
dtg=dat{:,5};
