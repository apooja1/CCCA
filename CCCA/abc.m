clc
clear all
% rng default
num=1000;
x=randn(num,1);
f1= @(x) x.^2; f2= @(y) y.^3;
y=f1(x);
z=f2(y);
% z=randn(num,1);

% Rho = [3 0 0 ;0 0.03 0;0 0 0.03];
% Z = lhsnorm([0 0 0], Rho, num);
% U = normcdf(Z,0,1);
% x=U(:,1); y=U(:,2); z=U(:,3);
rr=cos3d(x,y,z)
% y=f(x);
%  result= cos2d(x', y')

