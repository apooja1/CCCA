function [dKLMdB] = derKLB_BSS(B,X,alpha,copula)

% This function calculates the derivative of the KL diveregence with respect 
% to demixing matrix B

y=B*X';
y=y';
N=length(y(:,1)); 
p=length(y(1,:)); 
total=0;
[Fhat,h]=Fhatestimator(y);
[C_Yhat,H]=copulaestimator(Fhat);

[dc_dvu0,dc_du0,dc_dv0] = vder(Fhat,alpha,copula);
DC=[dc_du0  dc_dv0];
if copula=='G1'
C_alpha=copulapdf('Gumbel',[Fhat(:,1) Fhat(:,2)],alpha);
elseif copula=='G2'
C_alpha=copulapdf('Gaussian',[Fhat(:,1) Fhat(:,2)],alpha);
elseif  copula=='C'
C_alpha=copulapdf('Clayton',[Fhat(:,1) Fhat(:,2)],alpha);
elseif  copula=='F'
 C_alpha=copulapdf('Frank',[Fhat(:,1) Fhat(:,2)],alpha);
end

c=(1/sqrt((2*pi)));

%% when sources are independent

  for l=1:N              
t12(l,1)=(1/(N*h(1)))*sum((c*exp(-0.5*(((y(:,1)-y(l,1))./h(1)).^2))).*(X(:,2)-X(l,2)));
t21(l,1)=(1/(N*h(2)))*sum((c*exp(-0.5*(((y(:,2)-y(l,2))./h(2)).^2))).*(X(:,1)-X(l,1)));
t11(l,1)=(1/(N*h(1)))*sum((c*exp(-0.5*(((y(:,1)-y(l,1))./h(1)).^2))).*(X(:,1)-X(l,1)));
t22(l,1)=(1/(N*h(2)))*sum((c*exp(-0.5*(((y(:,2)-y(l,2))./h(2)).^2))).*(X(:,2)-X(l,2)));

  end
t01=[t12 t21];
t02=[t11 t22];
for n=1:N
t1=0; 
% 
 for ip=1:p
%      jp=1;
        for jp=1:p
                t3=(1/(N*h(ip)))*sum((c*exp(-0.5*(((y(:,ip)-y(n,ip))./h(ip)).^2))).*(X(:,jp)-X(n,jp)));

                t5=-2*(((c*exp(-0.5*(((Fhat(:,ip)-Fhat(n,ip))./H(ip)).^2)))).*((Fhat(:,ip)-Fhat(n,ip))./H(ip)));

                 if ip~=jp
                 t2=(t01(:,ip)-t3);
                 t4=(([c*exp(-0.5*(((Fhat(:,:)-Fhat(n,:))./H).^2))]));

                 t1=(1/(N*prod(H)))*sum((t4(:,1).*t5.*t4(:,2)).*t2*(1/H(ip)));
%                
%                  v1=(1/N)*(dc_dvu0(n,1))*dc_du0(n,1)*dc_dv0(n,1)*sum(t2.*t4(:,1).*t4(:,2));   %*
                 v1=(1/N)*(dc_dvu0(n,1))*prod(DC(n,:))*sum(t2.*t4(:,1).*t4(:,2));
%                    v1=(1/N)*(dc_dvu0(n,1))*sum(t2.*t4(:,1).*t4(:,2));
                 else 
                 t2=(t02(:,ip)-t3);   
                 t4=(([c*exp(-0.5*(((Fhat(:,:)-Fhat(n,:))./H).^2))]));

                 t1=(1/(N*prod(H)))*sum((t5).*t2*(1/H(ip)));
%                  v1=(1/N)*dc_du0(n,1)*dc_dv0(n,1)*sum(t2);
%                  v1=(1/N)*(dc_dvu0(n,1))*dc_du0(n,1)*dc_dv0(n,1)*sum(t2);  %* 
                 v1=(1/N)*(dc_dvu0(n,1))*prod(DC(n,:))*sum(t2);
%                  v1=(((1/N)*(dc_dvu0(n,1))*sum(t2)));
%                 v1= (1/N)*(dc_dvu0(n,1))*prod(DC(n,:))*sum(t2.*t4(:,1).*t4(:,2));
                 end

%                  t1=(1/(N*prod(H)));
                inde(ip,jp)=t1;
                inde0(ip,jp)=v1;
                t1=0;
        end
                 


                
      
%         jp=jp+1;

 end
         
   total=total+(inde./C_Yhat(n))-(inde0./C_alpha(n));
   inde=[];     

 end
  
 
  
dKLMdB=(1/N)*total;
end





%% when sources are dependent
