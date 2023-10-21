
clc 
clear all

n=60; p=10;

A=rand(n,p);

A(:,5)=2*A(:,4);  
A(:,2)=4*A(:,1);

[coeff, score]=pca(A);

Xcentered = score*coeff';
ndim=5;
[residuals,reconstructed] =pcares(zscore(A), ndim);

