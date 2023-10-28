function [KL] = KL_div(C_alpha,C_Yhat,N)

% KL=(1/N)* sum(log(C_Yhat./C_alpha));
KL=(1/N)*sum((log(C_Yhat)-log(C_alpha)));

end