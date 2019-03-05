%The second order ODE expressed as two coupled 1st order ODE.
function [xdot] = cfin_numerical(r,t)
    global l2;
    xdot=zeros(2,1);
    xdot(1)=t(2);
    xdot(2)=l2*t(1)-t(2)/r;
end