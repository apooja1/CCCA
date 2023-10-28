function [alpha] = Depedence_Param_estimation(COS,copula)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% inputs: 'C': Clayton, 
%         'F': Frank,
%         'G1': Gumbel, 
%         'G2': Gaussian, 

% load beta_Gau 
% load beta_Gum 
load beta_cla.mat
% load beta_fra.mat

h= [COS^2 COS 1];
% h= [COS^6 COS^5 COS^4 COS^3 COS^2 COS 1];

if copula=='C'
% h= [COS^6 COS^5 COS^4 COS^3 COS^2 COS 1];
alpha=h*beta_cla;

elseif copula=='F' 

alpha=h*beta_fra;

elseif copula=='G1' 

alpha=h*beta_Gum;

elseif copula=='G2' 
    
alpha=h*beta_Gau;
end


end

