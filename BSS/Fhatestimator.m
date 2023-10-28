function [Fhat,h] = Fhatestimator(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% h= bandwidth  (vector);
N=length(x(:,1));
p=length(x(1,:));




% pp = @(x1,m,s) (1/sqrt((2*pi)))*exp(-0.5*(((x1-m)./s).^2));   % probability density
%  pp = @(x1,m,s) exp(-((x1-m).^2)/(2*s.^2)) /(s*sqrt(2*pi));

% P_c

for j=1:p
for i=1:N
Mean = mean(x(i,j));
Standard_Deviation = std(x(:,j));
h(j)=(4/3)^(0.2)*N^(-0.2)*Standard_Deviation;
lims = [-Inf*ones(N,1), x(:,j)];
cp = normcdf(lims, Mean,sqrt(2*h(j)));
ee=cp(:,2) - cp(:,1);
Fhat(i,j)  = (1/N)*sum(ee);


end

end




end


