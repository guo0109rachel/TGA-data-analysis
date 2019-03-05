function C = NCQWO5()
%NCQWO5 Summary of this function goes here
%newton cotes quadrature weight
%for open boundary, n=5
%   Detailed explanation goes here
n=5;
C=zeros(1,n+1);
j=0:n;
%j(j==1)=[];
for i=0:n
    j(i+1)=[];
    %fun=@(s)(s-j)./(i-j));
    fun=@(s)(s-j(1)).*(s-j(2)).*(s-j(3)).*(s-j(4)).*(s-j(5))./(i-j(1))./(i-j(2))./(i-j(3))./(i-j(4))./(i-j(5));
    C(i+1)=1/7*integral(fun,0,n);
    j=0:n;
end



end

