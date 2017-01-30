function f = loglike_pk(y,part_vals,d)

% log-likelihood calculation for PK model
sigma2 = 0.1;
tau2 = 0.01;
D=400;
f = zeros(length(part_vals(:,1)),1);
theta = exp(part_vals);

for i = 1:length(d)
    
    mu = exp(-theta(:,1)*d(i)) - exp(-theta(:,2)*d(i));
    c = D./theta(:,3).*theta(:,2)./(theta(:,2)-theta(:,1));
    
    the_mean = c.*mu;
    v = 1+tau2/sigma2*mu.^2;
    the_var = sigma2*v;
    
    f = f - 0.5*log(2*pi*the_var) - 0.5./the_var.*(y(i) - the_mean).^2;
      
end

end