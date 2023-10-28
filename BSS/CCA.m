function [SNR,errorT1,B] = CCA(X,U0,copula,mode)
n=length(X(1,:));
p=length(X(:,1));  % enter the number of sources

maxIter=500;



%% Generate dependent samples from the FGM copula/ Clyaton/ Frank copula
% inputs: 'C': Clayton, 
%         'F': Frank,
%         'G1': Gumbel, 
%         'G2': Gaussian, 
% alpha=3;
% copula='G1';   
% U0 = copularnd('Gumbel',alpha,n);        
% s1=U0(:,1);s2=U0(:,2);

% 
%% generating independent sources
% [r1,r2] = RandStream.create('mrg32k3a','NumStreams',2);
% s1 = rand(r1,10000,1); 
% s2 = rand(r2,10000,1);
% corrcoef([s1,s2])
% COS=cosdv(s1,s2)
% s2=s2(1:n); s1=s1(1:n);
% U0=[s1 s2];
           

%% Forming the set of observations
% X=A*[s1'; s2']; %+0.001*randn(p,n);      

%% Intializations

iter=0;
error=1;
threshold=2;
mu1=0.01;
mu2=0.001;
% B0=-1*ones(p); B0(1,1)=1;  B0(2,2)=1; 
B0=eye(p);
SNRdiff=1;
% SNR=1;
y0=B0*X;

if copula=='G1'
     theta=5;
elseif copula=='G2'
     theta=0.5;
elseif ((copula=='C') || (copula=='F')) 
     theta=10;
end 

% Mdl = rica(y0',p,'Standardize',true);
% z = transform(Mdl,y0');
% z=z';
z=B0*X;

while ((error<threshold)&&~(SNRdiff<0))
    if mode=='D'
     theta1=theta-mu2*dKLMdtheta(z,theta,copula);  
    elseif mode=='I'
     theta1=1;   
    end
% if theta1<1
%     theta1=theta1+1;
% end
[dKLMdB,C_Yhat,C_alpha] = derKLB(z,X,(theta1),copula,mode);
dKLMdB;

B0;
B=B0+mu1*dKLMdB;
% B = permute(B,[2 1]);
B;
error1=(KL_div(C_alpha,C_Yhat,n));
error=real(ceil(log10(error1)-1))
B0=B;
theta=theta1;
iter=iter+1;
errorT1(iter,1)=error1;
thetaBig(iter,1)=theta;

z=B0*X;

if ~isnan(U0)
for i=1:p
SNR(iter,i)=10*log10(sum(U0(:,i).^2)/sum((z(i,:)-U0(:,i)').^2))
end
SNR;
end
if iter>maxIter
    break
end
if ~isnan(U0)
if iter>7
   SNRdiff=SNR(iter,1)-SNR(iter-1,1);
else 
SNRdiff=1;
end

end
end
end
% yhat=B*X;
% B
% Ainv=A^(-1) 
% figure()
% plot(errorT)
% title('Error')
% figure()
% plot(SNR')
% title('SNR')














