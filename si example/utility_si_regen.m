function u = utility_si_regen(d,J,prior_sims,alpha)


% SI model simulations
r = rand(J,2+50); % generate pseudo random numbers
theta(:,1) = exp(norminv(r(:,1),-3.5,sqrt(0.1024)));
theta(:,2) = exp(norminv(r(:,2),-4.5,sqrt(0.16)));

y=zeros(J,length(d));
for j = 1:J
    y(j,:) = simulate_SI(d,theta(j,1),theta(j,2),50,r(j,3:end));
end

u = zeros(J,1);

for j = 1:J
    % generate from prior predictive distribution
    r = rand(prior_sims,2+50); % generate pseudo random numbers
    theta_sim(:,1) = exp(norminv(r(:,1),-3.5,sqrt(0.1024)));
    theta_sim(:,2) = exp(norminv(r(:,2),-4.5,sqrt(0.16)));
    ysim=zeros(prior_sims,length(d));
    for i = 1:prior_sims
        ysim(i,:) = simulate_SI(d,theta_sim(i,1),theta_sim(i,2),50,r(i,3:end));
    end
    std_sim = std(ysim);
    std_sim_rep = repmat(std_sim,prior_sims,1);
    
    % ABC approximation
    rho = sum(abs(repmat(y(j,:),prior_sims,1) - ysim)./std_sim_rep,2);
    eps = quantile(rho,alpha);
    theta_post = theta_sim(rho<=eps,:);
    u(j) = -log(det(cov(theta_post)));
end

% approximate expected utility
u = mean(u);


end