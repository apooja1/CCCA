function [KLm] = KLdiveregence(B,X,alpha_est)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n=length(X(1,:)); p=length(X(:,1));
y=B*X;
[Fhat, Phat] = Fhatestimator(y');   % CDf estimation
u=Fhat(:,1);  v=Fhat(:,2);
C_Y = copulaestimator(Fhat);   % c_Yhat copula estimator  (non-parametric)
% C_Y=copulapdf('Frank',[u v],alpha_est);  %(parametric)
C_alpha=copulapdf('Frank',[u v],alpha_est);  
KLm=(1/n)*sum(log(C_Y./C_alpha));   %when sources are dependent







% KLm=-(1/n)*sum(log(C_Y));   %when sources are independent
end

