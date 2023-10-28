function [output] = dKLMdtheta(y,theta,copula)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% y=B*X;
y=y';
N=length(y(:,1)); 
% p=length(y(1,:)); 
alpha=theta;
[Fhat,~]=Fhatestimator(y);
% [C_Yhat,H]=copulaestimator(Fhat);
 
if copula=='G1'        % Gumbel copula
% syms u v alpha 
% f=exp(-(((-log(u)).^(alpha)+((-log(v)).^(alpha)))).^(1/alpha)); %differentiation of log of f
% dc=diff(diff(f,u),v);
% diff(dc,alpha)
C_alpha=copulapdf('Gumbel',[Fhat(:,1) Fhat(:,2)],abs(alpha));
dc_dalpha= @(u,v,alpha) (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*log(-log(u))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(2/alpha - 2))/(u*v) - (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((2*log((-log(u))^alpha + (-log(v))^alpha)*((-log(u))^alpha + (-log(v))^alpha)^(2/alpha - 2))/alpha^2 - ((-log(u))^alpha + (-log(v))^alpha)^(2/alpha - 3)*(log(-log(u))*(-log(u))^alpha + log(-log(v))*(-log(v))^alpha)*(2/alpha - 2)))/(u*v) + (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*log(-log(v))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(2/alpha - 2))/(u*v) - (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*(1/alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/(u*v) + (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((log((-log(u))^alpha + (-log(v))^alpha)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))/alpha^2 - (((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 1)*(log(-log(u))*(-log(u))^alpha + log(-log(v))*(-log(v))^alpha))/alpha)*((-log(u))^alpha + (-log(v))^alpha)^(2/alpha - 2))/(u*v) + (exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/(alpha*u*v) + (alpha*exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*((log((-log(u))^alpha + (-log(v))^alpha)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/alpha^2 - (1/alpha - 2)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 3)*(log(-log(u))*(-log(u))^alpha + log(-log(v))*(-log(v))^alpha))*(1/alpha - 1))/(u*v) - (alpha*exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*log(-log(u))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*(1/alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/(u*v) - (alpha*exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*log(-log(v))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*(1/alpha - 1)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/(u*v) - (alpha*exp(-((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))*(-log(u))^(alpha - 1)*(-log(v))^(alpha - 1)*(1/alpha - 1)*((log((-log(u))^alpha + (-log(v))^alpha)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha))/alpha^2 - (((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 1)*(log(-log(u))*(-log(u))^alpha + log(-log(v))*(-log(v))^alpha))/alpha)*((-log(u))^alpha + (-log(v))^alpha)^(1/alpha - 2))/(u*v); 

elseif copula=='G2'        % Gaussian copula
C_alpha=copulapdf('Gaussian',[Fhat(:,1) Fhat(:,2)],alpha);
dc_dalpha= @(u,v,alpha) (8*alpha^2*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2)))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)^2) - (2*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2)))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)) - (2*alpha^2*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2)))/((1 - alpha^2)^(3/2)*(2*alpha^2 - 2)) + (2*alpha*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*((2*u*v - 2*alpha*(u^2 + v^2))/(2*alpha^2 - 2) + (4*alpha*(alpha^2*(u^2 + v^2) - 2*alpha*u*v))/(2*alpha^2 - 2)^2))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)) + (exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*(- 2*v*alpha^2 + 2*u*alpha)*(2*v - 4*alpha*u))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)^2) + (exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*(- 2*u*alpha^2 + 2*v*alpha)*(2*u - 4*alpha*v))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)^2) - (exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*(- 2*v*alpha^2 + 2*u*alpha)*(- 2*u*alpha^2 + 2*v*alpha)*((2*u*v - 2*alpha*(u^2 + v^2))/(2*alpha^2 - 2) + (4*alpha*(alpha^2*(u^2 + v^2) - 2*alpha*u*v))/(2*alpha^2 - 2)^2))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)^2) - (8*alpha*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*(- 2*v*alpha^2 + 2*u*alpha)*(- 2*u*alpha^2 + 2*v*alpha))/((1 - alpha^2)^(1/2)*(2*alpha^2 - 2)^3) + (alpha*exp((alpha^2*(u^2 + v^2) - 2*alpha*u*v)/(2*alpha^2 - 2))*(- 2*v*alpha^2 + 2*u*alpha)*(- 2*u*alpha^2 + 2*v*alpha))/((1 - alpha^2)^(3/2)*(2*alpha^2 - 2)^2);


elseif copula=='C'    %Clayton
%    syms u v alpha1
%    f=log(((u)^(-alpha1)+(v)^(-alpha1)-1)^(-1/alpha1)); 
%    dc=diff(diff(f,u), v); 
%    diff(dc,alpha1)
C_alpha=copulapdf('Clayton',[Fhat(:,1) Fhat(:,2)],alpha);
% dc_dalpha= @(u,v,alpha1) (1/alpha1 + 1)/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) - 1/(alpha1*u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) + (alpha1*(1/alpha1 + 1)*(log(1/u^alpha1 + 1/v^alpha1 - 1)/(alpha1^2*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) + ((log(u)/u^alpha1 + log(v)/v^alpha1)*(1/alpha1 + 2))/(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 3)))/(u^(alpha1 + 1)*v^(alpha1 + 1)) - (alpha1*log(u)*(1/alpha1 + 1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) - (alpha1*log(v)*(1/alpha1 + 1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2));
  dc_dalpha= @(u,v,alpha1) ((log(1/u^alpha1 + 1/v^alpha1 - 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 1))/alpha1^2 + (log(u)/u^alpha1 + log(v)/v^alpha1)*(1/alpha1 - 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 2))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 1)) - ((log(1/u^alpha1 + 1/v^alpha1 - 1)/(alpha1^2*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 1)) + ((log(u)/u^alpha1 + log(v)/v^alpha1)*(1/alpha1 + 1))/(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2))*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 1))/(u^(alpha1 + 1)*v^(alpha1 + 1)) - (1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1)/(alpha1*u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) + ((1/alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) + (log(u)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 1)) + (log(v)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 1)) + (alpha1*(1/alpha1 + 1)*(log(1/u^alpha1 + 1/v^alpha1 - 1)/(alpha1^2*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) + ((log(u)/u^alpha1 + log(v)/v^alpha1)*(1/alpha1 + 2))/(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 3))*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1))/(u^(alpha1 + 1)*v^(alpha1 + 1)) - (alpha1*(1/alpha1 + 1)*(((log(u)/u^alpha1 + log(v)/v^alpha1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 - 1))/alpha1 + (log(1/u^alpha1 + 1/v^alpha1 - 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1))/alpha1^2))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) - (alpha1*log(u)*(1/alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2)) - (alpha1*log(v)*(1/alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1))/(u^(alpha1 + 1)*v^(alpha1 + 1)*(1/u^alpha1 + 1/v^alpha1 - 1)^(1/alpha1 + 2));

elseif copula=='F'   %Frank   
%     syms u v alpha1
%     f=(-1/alpha1)*log(1+((exp(-alpha1*u-1)*exp(-alpha1*v-1))/(exp(-alpha1)-1)));
%     dc=diff(diff(f,u), v); 
%     diff(dc,alpha1)
C_alpha=copulapdf('Frank',[Fhat(:,1) Fhat(:,2)],alpha);
dc_dalpha= @(u,v,alpha1) (exp(- 2*alpha1*u - 2)*exp(- 2*alpha1*v - 2))/((exp(-alpha1) - 1)^2*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^2) - (exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/((exp(-alpha1) - 1)*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)) + (alpha1*u*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/((exp(-alpha1) - 1)*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)) - (2*alpha1*u*exp(- 2*alpha1*u - 2)*exp(- 2*alpha1*v - 2))/((exp(-alpha1) - 1)^2*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^2) + (alpha1*v*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/((exp(-alpha1) - 1)*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)) - (2*alpha1*v*exp(- 2*alpha1*u - 2)*exp(- 2*alpha1*v - 2))/((exp(-alpha1) - 1)^2*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^2) - (alpha1*exp(- alpha1*u - 1)*exp(- alpha1*v - 1)*((u*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) - (exp(-alpha1)*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1)^2 + (v*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1)))/((exp(-alpha1) - 1)*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^2) + (2*alpha1*exp(- 2*alpha1*u - 2)*exp(- 2*alpha1*v - 2)*((u*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) - (exp(-alpha1)*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1)^2 + (v*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1)))/((exp(-alpha1) - 1)^2*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^3) - (alpha1*exp(-alpha1)*exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/((exp(-alpha1) - 1)^2*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)) + (2*alpha1*exp(-alpha1)*exp(- 2*alpha1*u - 2)*exp(- 2*alpha1*v - 2))/((exp(-alpha1) - 1)^3*((exp(- alpha1*u - 1)*exp(- alpha1*v - 1))/(exp(-alpha1) - 1) + 1)^2);
end

for n=1:N
 Out(n,1)= dc_dalpha(Fhat(n,1),Fhat(n,2),alpha);
end
%  output=-(1/N)*sum(Out./C_alpha);
 output=-(1/N)*sum(Out);
end

