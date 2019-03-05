%% This file is only appropriate for TGA_redox cycles activation calculation
% Only be available when iron oxide as active materials and oxygen carrier
% Change the conversion calculation when used as other materials
% The 100% conversion means Fe3O4
% afterward change the output line
% Third, check cycles,thres,flush time and oxidation time of each
% experiments and revise the parameters section below

clear
clc
close('all')


% Parameters
n_cycle=input('how many cycles you have for your reaction?...'); % number of redox cycles

%Temprature for the reaction
Temp=input('What is your temperature of reaction?...use ";" to separate different temperature...');
T=1./(Temp+273.15);
n_temp=length(Temp);

%Construct the data management
data=zeros(n_cycle,19,n_temp);

%data read from the excel
[FileName,PathName,~]=uigetfile('*.xls*','select the data for activation energy calculation');
for i=1:n_temp
    data(:,:,i)=xlsread([PathName,FileName],num2str(Temp(i)));
end

%Calculate the m average value
for j=1:n_cycle
    m_average(:,j)=sum(data(j,2,:))/n_temp;
end

%polyfit to calculate m and r2
for p=1:n_cycle
    for q=1:8
        for o=1:n_temp
        k(o)=data(p,(q+1),o);
        end
        [b,~,~,~,s]=regress(log(k'),[ones(n_temp,1),T']);
        E(p,q)=b(2)*(-8.314);
        ER2(p,q)=s(1);
    end
end
%when m=0.54 get E(1) for three-dimensional diffusion

%when m=0.57 get E(2) for two-dimensional diffusion

%when m=0.62 get E(3) for one-dimensional diffusion

%when m=1 get E(4) for first order reaction;

%when m=1.07 get E(5) for sphere Phase-boundary Controlled

%when m=1.11 get E(6) for cylinder Phase-boundary Controlled

%when m=2 get E(7) for two-dimensional growth of nuclei

%when m=3 get E(8) for three-dimensional growth of nuclei










