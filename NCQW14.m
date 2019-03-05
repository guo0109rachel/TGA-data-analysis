function C = NCQW14()
%NCQW14 quadrature weights for n=14
%   Detailed explanation goes here
n=14;
C=zeros(1,n+1);
j=0:n;
%j(j==1)=[];
for i=0:n
    j(i+1)=[];
    %fun=@(s)(s-j)./(i-j));
    fun=@(s)(s-j(1)).*(s-j(2)).*(s-j(3)).*(s-j(4)).*(s-j(5)).*(s-j(6)).*(s-j(7)).*(s-j(8)).*(s-j(9)).*(s-j(10)).*(s-j(11)).*(s-j(12)).*(s-j(13)).*(s-j(14))./(i-j(1))./(i-j(2))./(i-j(3))./(i-j(4))./(i-j(5))./(i-j(6))./(i-j(7))./(i-j(8))./(i-j(9))./(i-j(10))./(i-j(11))./(i-j(12))./(i-j(13))./(i-j(14));
    C(i+1)=1/7*integral(fun,0,n);
    j=0:n;
end

end

