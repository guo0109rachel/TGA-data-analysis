function [Weightper,Conversion] = TGA_convcalc_withOindopant(rawdata,dopant,percent)
%percent means the dopant percentage
%weight per is the weight percentage of total weight
%Conversion is based on oxygen conversion
Maxweight=max(rawdata);
Weightper=rawdata/Maxweight*100;
%input the relative weight of atom
O=15.999; Fe=55.845; La=138.905; Co=58.993; Ni=58.693; Cu=63.546;
%dopant depends on different dopant 0=non,1=La,2=Co,3=Ni,4=Cu
%dopant=0 means no dopant, only Fe2O3
if dopant==0
    Conversion=(100-Weightper)/((O*3/(O*3+Fe*2))*100)*100;
    %dopant=1 means La dopant
else if dopant==1
        WeightFe2O3=100*(2*Fe+3*O)/(percent/100*(2*La+3*O)+(O*3+Fe*2));
        WeightLa2O3=100-WeightFe2O3;
        Conversion=(100-Weightper)/((O*3/(O*3+Fe*2))*WeightFe2O3+(O*3/(O*3+La*2))*WeightLa2O3)*100;
    %dopant=2 means Co dopant
    else if dopant==2
        WeightFe2O3=100*(2*Fe+3*O)/(percent/100*(2*Co+2*O)+(O*3+Fe*2));
        WeightCoO=100-WeightFe2O3;
        Conversion=(100-Weightper)/((O*3/(O*3+Fe*2))*WeightFe2O3+(O/(O+Co))*WeightCoO)*100;
    %dopant=3 means Ni dopant
        else if dopant==3
        WeightFe2O3=100*(2*Fe+3*O)/(percent/100*(2*Ni+2*O)+(O*3+Fe*2));
        WeightNiO=100-WeightFe2O3;
        Conversion=(100-Weightper)/((O*3/(O*3+Fe*2))*WeightFe2O3+(O/(O+Ni))*WeightNiO)*100;
    %dopant=4 means Cu dopant
            else if dopant==4
        WeightFe2O3=100*(2*Fe+3*O)/(percent/100*(2*Cu+2*O)+(O*3+Fe*2));
        WeightCuO=100-WeightFe2O3;
        Conversion=(100-Weightper)/((O*3/(O*3+Fe*2))*WeightFe2O3+(O/(O+Cu))*WeightCuO)*100;
                end
            end
        end
    end
end
end