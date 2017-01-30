function u = utility_discrim(d,J,m_sim,prior_sims,ysim,std_sim_rep,alpha)

% death model simulations
r1 = rand(J,1+50); % generate pseudo random numbers
theta1 = exp(norminv(r1(:,1),-0.48,0.15));

y1=zeros(J,length(d));
for j = 1:J
    y1(j,:) = simulate_death(d,theta1(j,1),50,r1(j,2:end));
end

% SI model simulations
r2 = rand(J,2+50); % generate pseudo random numbers
theta2(:,1) = exp(norminv(r2(:,1),-1.1,0.2));
theta2(:,2) = exp(norminv(r2(:,2),-4.5,0.6));

y2=zeros(J,length(d));
for j = 1:J
    y2(j,:) = simulate_SI(d,theta2(j,1),theta2(j,2),50,r2(j,3:end));
end

% ABC approximations
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