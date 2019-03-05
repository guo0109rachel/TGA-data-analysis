%This file is only used combine with *TGATPR_MS_write program
function TGA_MSresults = TGA_MSfunc_direct(TGAdata,TMS,MSgasraw,dopant,percent,Helium_flow,delay)
%Handle the data for TGA_MS results for TPR reaction

%Find the highest temperature in TGA
[~,rows]=max(TGAdata(:,3));

%Delete the useless data
TGAcutdata=TGAdata(1:rows,:);

%Find the whole time
tTGA=TGAcutdata(rows,2)*60-delay;%unit is second
tMS=max(TMS);%unit is second
if tTGA>tMS
    tTOTAL=floor(tMS);
else tTOTAL=floor(tTGA);
end

%Obtain some useful raw TGA data
TGA_rawtime=TGAcutdata(:,2)*60;%unit is second
TGA_rawtemp=TGAcutdata(:,3);
TGA_rawweight=TGAcutdata(:,4);
TGA_rawDTG=TGAcutdata(:,5);
Maxweight=max(TGA_rawweight);
TGA_rawweightper=TGA_rawweight/Maxweight*100;

%set the timeline
t=0:0.5:tTOTAL;
t=t';

%Do the line fitting for TGA
TGA_temp=spline(TGA_rawtime,TGA_rawtemp,t);
TGA_weight=spline(TGA_rawtime,TGA_rawweight,t);
TGA_DTG=spline(TGA_rawtime,TGA_rawDTG,t);
[TGA_weightper,Conversion]=TGA_convcalc(TGA_weight,dopant,percent);

%CO2 is the first one, Helium is the second, CO and H2 and then CH4;
%Obtain the gas data
CO2raw=MSgasraw(:,1);
Heliumraw=MSgasraw(:,2);
COraw=MSgasraw(:,3);
H2raw=MSgasraw(:,4);
CH4raw=MSgasraw(:,5);

%Do the line fitting for MS
CH4=spline(TMS,CH4raw,t);
CO2=spline(TMS,CO2raw,t);
CO=spline(TMS,COraw,t);
H2=spline(TMS,H2raw,t);
Helium=spline(TMS,Heliumraw,t);

%Do the flow rate calculation
Totalflow=Helium_flow./Helium*100;
CH4_flow=Totalflow.*CH4/100;
CO2_flow=Totalflow.*CO2/100;
CO_flow=Totalflow.*CO/100;
H2_flow=Totalflow.*H2/100;

%Do the concentration without Helium
FlowwithoutHe=CH4_flow+CO2_flow+CO_flow+H2_flow;
CH4_percent=CH4_flow./FlowwithoutHe;
CO2_percent=CO2_flow./FlowwithoutHe;
CO_percent=CO_flow./FlowwithoutHe;
H2_percent=H2_flow./FlowwithoutHe;

%assemble the final data
TGA_MSresults=table(t,t/60,TGA_temp,TGA_weight,TGA_weightper,TGA_DTG,Conversion,CO2,Helium,CO,H2,CH4,...
    CO2_flow,CO_flow,H2_flow,CH4_flow,CO2_percent,CO_percent,H2_percent,CH4_percent);

end

