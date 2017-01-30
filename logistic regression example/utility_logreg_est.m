function u = utility_logreg_est(X,J,sims,mu_prior,std_prior)

u = zeros(J,1);
[n,~] = size(X);

a = repmat(mu_prior,sims,1);
b = repmat(std_prior,sims,1);

for j = 1:J
    
    % draws from prior distribution for importance sampling
    beta_prior = normrnd(a,b);
    
    % draw potential dataset from prior predictive
    beta = normrnd(mu_prior,std_prior);
    p = 1./(1 + exp(-beta(1) - beta(2)*X(:,1) - beta(3)*X(:,2) - beta(4)*X(:,3) - beta(5)*X(:,4)));
    y = zeros(n,1);
    y(rand(n,1)<p) = 1;
    
    % importance sampling approximation of posterior and utility u(d,y)
    f = loglike_logreg(y,beta_prior,X);
    w = f - max(f);
    w = exp(w);
    w = w/sum(w);
    u(j) = sum(w(w>0).*f(w>0)) + log(sims) - logsumexp(f);
    
end

% approximate the expected utility
u = mean(u);

end


