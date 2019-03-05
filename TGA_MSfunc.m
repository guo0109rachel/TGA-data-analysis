function TGA_MSresults = TGA_MSfunc(TGAdata,TMS,MSgasraw,dopant,percent,delay)
%Handle the data for TGA_MS results for TPR reaction

%Find the highest temperature in TGA
[max_temp,rows]=max(TGAdata(:,3));

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

%Obtain the gas data
CO2raw=MSgasraw(:,1);
Heliumraw=MSgasraw(:,4);
COraw=MSgasraw(:,7);
H2raw=MSgasraw(:,10);
CH4raw=MSgasraw(:,13);

%Do the line fitting for MS
CH4=spline(TMS,CH4raw,t);
CO2=spline(TMS,CO2raw,t);
CO=spline(TMS,COraw,t);
H2=spline(TMS,H2raw,t);
Helium=spline(TMS,Heliumraw,t);

%assemble the final data
TGA_MSresults=[t,t/60,TGA_temp,TGA_weight,TGA_weightper,TGA_DTG,Conversion,CO2,Helium,CO,H2,CH4];

end

