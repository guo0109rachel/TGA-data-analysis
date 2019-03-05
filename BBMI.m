function p = BBMI( x,y )
%BBMI Summary of this function goes here
%the monochromatic intensity of the black body
%x represents the wavelength, lambda, meter
%y represents the temperature, T, kelvin

C1=1.19e-16;%unit: W/m^2
C2=1.4388e-2;%unit: m*K

p=C1/x^5./(exp(C2/x./y)-1);

end

