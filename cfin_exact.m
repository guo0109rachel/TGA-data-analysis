% The exact solution in terms of Bessel functions.
function [theta]=cfin_exact(r)
global lambda;
global R1;
global R2;

l1 = lambda * R1;
l2 = lambda * R2;

s=besseli(0,l1)*besselk(1,l2)+besseli(1,l2)*besselk(0,l1);
c=besselk(1,l2)/s;
d=besseli(1,l2)/s;

lr = lambda * r;
theta = c * besseli(0,lr)+d*besselk(0,lr);
end
