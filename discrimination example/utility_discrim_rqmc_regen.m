function u = utility_discrim_rqmc_regen(d,J,prior_sims,alpha)

% generate from prior predictive of death model
r1 = rand(prior_sims,1+50); % generate pseudo random numbers
theta_sim1 = exp(norminv(r1(:,1),-0.48,0.15));
ysim1=zeros(prior_sims,length(d));
for j = 1:prior_sims
    ysim1(j,:) = simulate_death(d,theta_sim1(j,1),50,r1(j,2:end));
end

% generate from prior predictive of SI model
r2 = rand(prior_sims,2+50); % generate pseudo random numbers
theta_sim2(:,1) = exp(norminv(r2(:,1),-1.1,0.2));
theta_sim2(:,2) = exp(norminv(r2(:,2),-4.5,0.6));
ysim2=zeros(prior_sims,length(d));
for j = 1:prior_sims
    ysim2(j,:) = simulate_SI(d,theta_sim2(j,1),theta_sim2(j,2),50,r2(j,3:end));
end

ysim = [ysim1; ysim2];
m_sim = [ones(prior_sims,1); 2*ones(prior_sims,1)];

std_sim = std(ysim);
std_sim_rep = repmat(std_sim,2*prior_sims,1);

% death model simulations
r1 = gen_Sobol(ceil(log2(J)),1+50)'; % generate randomised QMC numbers
r1 = r1(1:J,:);
theta1 = exp(norminv(r1(:,1),-0.48,0.15));

y1=zeros(J,length(d));
for j = 1:J
    y1(j,:) = simulate_death(d,theta1(j,1),50,r1(j,2:end));
end

% SI model simulations
r2 = gen_Sobol(ceil(log2(J)),2+50)'; % generate randomised QMC numbers
r2 = r2(1:J,:);
theta2(:,1) = exp(norminv(r2(:,1),-1.1,0.2));
theta2(:,2) = exp(norminv(r2(:,2),-4.5,0.6));

y2=zeros(J,length(d));
for j = 1:J
    y2(j,:) = simulate_SI(d,theta2(j,1),theta2(j,2),50,r2(j,3:end));
end

u1 = zeros(J,1);

for j = 1:J
    rho = sum(abs(repmat(y1(j,:),2*prior_sims,1) - ysim)./std_sim_rep,2);
    eps = quantile(rho,alpha);
    m_post = m_sim(rho<=eps);
    u1(j) = log(mean(m_post == 1));
end

u2 = zeros(J,1);

for j = 1:J
    rho = sum(abs(repmat(y2(j,:),2*prior_sims,1) - ysim)./std_sim_rep,2);
    eps = quantile(rho,alpha);
    m_post = m_sim(rho<=eps);
    u2(j) = log(mean(m_post == 2)); 
end

% estimate expected utility
u = 0.5*u1 + 0.5*u2;

u = mean(u);

end