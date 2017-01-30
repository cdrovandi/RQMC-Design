function f = loglike_logreg(y,beta,X)

%beta = (prior_lhs + exp(beta_t).*prior_rhs)./(1 + exp(beta_t));

[N,~] = size(beta);


f = zeros(N,1);

for i = 1:length(y)
    p = 1./(1 + exp(-beta(:,1) - beta(:,2)*X(i,1) - beta(:,3)*X(i,2) - beta(:,4)*X(i,3) - beta(:,5)*X(i,4)));
    if y(i) == 1
        f = f + log(p);
    else
        f = f + log(1-p);
    end
end


end