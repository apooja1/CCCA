clc 
clear all 

close all
% format long
% rng default

p=2;  % enter the number of sources
maxIter=100;

A=[1 0.4; 0.4 1];

run('data_4mc_2area.m');
d1=output_data(:,2); 
d2=output_data(:,3); 
d3=output_data(:,4); 
d4=output_data(:,5); 

n=length(d1);

%% Generate dependent samples from the FGM copula/ Clyaton/ Frank copula
% inputs: 'C': Clayton,  % 0.1:1
%         'F': Frank,    % 0.1:1
%         'G1': Gumbel,  %1:2
%         'G2': Gaussian,%0.09 to 0.99 
       
s1=d3; % update source signals as rotor speed 1  
s2=d4; % source 2
mode='D';

%% generating independent sources
% [r1,r2] = RandStream.create('mrg32k3a','NumStreams',2);
% s1 = rand(r1,10000,1); 
% s2 = rand(r2,10000,1);
% corrcoef([s1,s2])
% COS=cosdv(s1,s2)
% s2=s2(1:n); s1=s1(1:n);
U0=[s1 s2];
           
copula='C'; 

%% Forming the set of observations
X=A*[s1'; s2']+0.001*randn(p,n);      
COS=cosdv(X(1,:)',X(2,:)');
%% Intializations
iter=0;
error=1;
threshold=2;
mu1=0.01; %stepsize for B
mu2=0.001; %stepsize for theta
B0=eye(p);
SNRdiff=1;
y0=B0*X;
z=B0*X;
theta1=0.7;
% Mdl = rica(y0',p);
% z = transform(Mdl,y0');
% z=z';
while ((error<threshold)&&~(SNRdiff<0))
% Mdl = rica(z',p,"Standardize",true);
% z = transform(Mdl,z');
% z=(z-min(z))./(max(z)-min(z));
% z=z';
theta1=cosdv(z(1,:)',z(2,:)');  % feed nX1 vectors
theta1=Depedence_Param_estimation(theta1,copula)
theta1
% theta1=theta1-mu2*dKLMdtheta(z,abs(theta1),copula);  % feed pXn z matrix
% (we dont need this step)

[dKLMdB,C_Yhat,C_alpha] = derKLB(z,X,abs(theta1),copula,mode);
dKLMdB
error1=(KL_div(C_alpha,C_Yhat,n));
error=real(ceil(log10(error1)-1))
B=B0-mu1*dKLMdB;
B
% error=norm(B-B0)
B0=B;
iter=iter+1
errorT(iter,1)=error1;
yhat=B*X;
y0=B*X;
z=B*X;

for i=1:p
SNR(iter,i)=10*log10(sum(U0(:,i).^2)/sum((z(i,:)-U0(:,i)').^2));
end
SNR
if iter>maxIter
    break
end

if iter>2
    SNRdiff=SNR(iter,1)-SNR(iter-1,1);
else 
SNRdiff=1;
end

end
B_CCCA=B;
%%
[SNR_CCA,errorT_CCA,B_CCA] = CCA(X,U0,copula,mode);

%%
% yhat=B*X;
% B

Ainv=A^(-1) 
B_CCCA
B_CCA
figure()
plot(errorT,'linewidth',2)
hold on 
plot(errorT_CCA,'linewidth',2)
title('Error')
xlabel('Iterations')
ylabel('Error')
legend('CCCA','CCA')
set(gcf,'color','w')


SNR(end,:)=[]; SNR_CCA(end,:)=[]; 
figure()
plot(SNR,'linewidth',2)
hold on 
plot(SNR_CCA,'linewidth',2,'linestyle','-.')
hold off
xlabel('Iterations')
ylabel('SNR')
legend('CCCA (y_{1})','CCCA (y_{2})','CCA (y_{1})','CCA (y_{2})')
set(gcf,'color','w')
set(gca,'FontSize',10,'fontweight','bold')
set(gca, 'LineWidth', 1.5)
legend boxoff 
% Adjust font
set(gca, 'FontName', 'Helvetica')
