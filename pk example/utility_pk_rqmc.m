function f = utility_pk_rqmc(d,J)

sigma2 = 0.1;
tau2 = 0.01;
D=400;
u = zeros(J,1);
r = gen_Sobol(ceil(log2(J)),3+length(d))'; % generate randomised QMC numbers

for j = 1:J
    % simulate from prior predictive based on RQMC numbers
	prior = norminv(r(j,1:3),[log(0.1) log(1) log(20)], sqrt(0.05));
    theta = exp(prior);
    
    mu = exp(-theta(1)*d) - exp(-theta(2)*d);
    c = D/theta(3)*theta(2)/(theta(2)-theta(1));
    
    the_mean = c*mu;
    v = 1+tau2/sigma2*mu.^2;
    the_var = sigma2*v;

	y = norminv(r(j,4:end),the_mean,sqrt(the_var));
    
    % compute laplace approximation and approximate utility u(d,y)
    u(j) = laplace_pk(y,d,prior);
    
end

% approximate expected utility
f = mean(u);

end