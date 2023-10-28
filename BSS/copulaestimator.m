function [K,H] = copulaestimator(Fhat)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here  
% u:independent vector 

K=[];
N=length(Fhat(:,1));
p=length(Fhat(1,:));
c=(1/sqrt((2*pi)));
for i=1:p
H(i)=(4/(p+2))^(1/(p+4))*std(Fhat(:,i));
end




for m=1:N
  e1= (c*exp(-0.5*(((Fhat(:,:)-Fhat(m,:))./H).^2)));
 
  K(m,1) =(1/(N*prod(H)))*sum(prod(e1,2));
    
end
    
    
    

end

