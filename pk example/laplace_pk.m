function kld =  laplace_pk(y,d,theta)

%% SMC with constraints code

[x,~,~,~,~,hessian] = fminunc(@(x)-loglike_pk(y,x,d) - sum(log(normpdf(x,[log(0.1) log(1) log(20)],sqrt(0.05)))),theta,optimset('Display','off'));

cov_prior = 0.05*eye(3);
cov_post = inv(hessian);

mu_prior = [log(0.1) log(1) log(20)];
mu_post = x;

kld = 0.5*(trace(inv(cov_prior)*cov_post) + (mu_post - mu_prior)*inv(cov_prior)*(mu_post - mu_prior)' - 3 + log(det(cov_prior)) - log(det(cov_post)));

end


