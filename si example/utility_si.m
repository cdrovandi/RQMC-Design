function u = utility_si(d,J,theta_sim,prior_sims,ysim,std_sim_rep,alpha)

% SI model simulations
r = rand(J,2+50); % generate pseudo random numbers
theta(:,1) = exp(norminv(r(:,1),-3.5,sqrt(0.1024)));
theta(:,2) = exp(norminv(r(:,2),-4.5,sqrt(0.16)));

y=zeros(J,length(d));
for j = 1:J
    y(j,:) = simulate_SI(d,theta(j,1),theta(j,2),50,r(j,3:end));
end

% ABC approximations
u = zeros(J,1);

for j = 1:J
    rho = sum(abs(repmat(y(j,:),prior_sims,1) - ysim)./std_sim_rep,2);
    eps = quantile(rho,alpha);
    theta_post = theta_sim(rho<=eps,:);
    u(j) = -log(det(cov(theta_post)));
end

u = mean(u);


end