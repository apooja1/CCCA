clc 
clear all 
close all
% format long
rng default
n=500;
p=2;  % enter the number of sources

maxIter=500;
A=[1 0.8 ; 0.8 1];     %mixing matrix


%% Generate dependent samples from the FGM copula/ Clyaton/ Frank copula
% inputs: 'C': Clayton,  % 0.1:1
%         'F': Frank,    % 0.1:1
%         'G1': Gumbel,  %1:2
%         'G2': Gaussian,%0.09 to 0.99 
alpha=10;
copula='C';   
U0 = copularnd('Clayton',alpha,n);        
s1=U0(:,1);s2=U0(:,2);
mode='D';


%% Forming the set of observations
X=A*[s1'; s2']; %+0.001*randn(p,n);      
COS=cosdv(X(1,:)',X(2,:)');
[alpha_est] = Depedence_Param_estimation(COS,copula); 
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

if copula=='G1'
     theta1=20;
elseif copula=='G2'
     theta1=0.5;
elseif ((copula=='C') || (copula=='F')) 
     theta1=10;
end


while ((error<threshold)&&~(SNRdiff<0))
theta1=cosdv(z(1,:)',z(2,:)');  % feed nX1 vectors
theta1=Depedence_Param_estimation(theta1,copula)
theta1
% theta1=theta1-mu2*dKLMdtheta(z,theta1,copula);  % feed pXn z matrix
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
%%
[SNR_CCA,errorT_CCA,B_CCA] = CCA(X,U0,copula,mode);

%%
% yhat=B*X;
% B
B_CCA
Ainv=A^(-1) 
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





