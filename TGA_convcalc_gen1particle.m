function [Weightper,Conversion,index] = TGA_convcalc_gen1particle(rawdata,Support,initial_weight,high_T_weight)
%percent means the dopant percentage
%weight per is the weight percentage of total weight
%Conversion is based on oxygen conversion
[Maxweight,index]=max(rawdata);
Full_Weight=(Maxweight-(high_T_weight-initial_weight))*(100-Support)/100;
Weightper = (Maxweight-rawdata)/Full_Weight*100;
%input the relative weight of atom
O=15.999; Fe=55.845;
    Conversion=Weightper/((O*3/(O*3+Fe*2))*100)*100;
end