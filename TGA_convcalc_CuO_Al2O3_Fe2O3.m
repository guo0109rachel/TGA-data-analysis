%This function is only based on Weight percentage calculation!!!!
function [Weightper,Conversion] = TGA_convcalc_CuO_Al2O3_Fe2O3(rawdata,Fe2O3,CuO)
%Fe2O3 means the wt% of Fe2O3
%CuO means the wt% of CuO
%Al2O3 is set to 1-CuO-Fe2O3
%Conversion is based on oxygen conversion
Maxweight=max(rawdata);
Weightper=rawdata/Maxweight*100;
%input the relative weight of atom
O=15.999; Fe=55.845; Cu=63.546;
%conversion calculation
Conversion=(100-Weightper)/(Fe2O3*(3*O/(3*O+2*Fe))+CuO*(O/(Cu+O)))*100;
end