
%% Compiling of C code (requires that Matlab has access to a suitable C compiler)
mex simulate_death.c;
mex simulate_SI.c;

%% DISCRIMINATION

prior_sims = 100000; % number of ABC importance samples
alpha = 0.01; % proportion of samples to "keep" in ABC approximation 
d = [5 10]; % test design

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

J=[5 10 20 30 40 50 60 70 80 90 100]; % number of Monte Carlo draws to test for
repeats = 100;

% loops can be performed using parfor (parallel computing) if desired
for j = 1:length(J)
    J(j)
	
	% MC Integration with standard pseudo random numbers
    u = zeros(repeats,1);
    for r = 1:repeats
        u(r) = utility_discrim(d,J(j),m_sim,prior_sims,ysim,std_sim_rep,alpha);
    end
    u_std(j) = std(u);
    u_mean(j) = mean(u);
    
	% MC Integration with randomised quasi Monte Carlo random numbers
    for r = 1:repeats
        u(r) = utility_discrim_rqmc(d,J(j),m_sim,prior_sims,ysim,std_sim_rep,alpha);
    end
    u_stdr(j) = std(u);
    u_meanr(j) = mean(u);
    
end

% plot results
plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')

%% DISCRIMINATION regeneration

prior_sims = 100000; % number of ABC importance samples
alpha = 0.01; % proportion of samples to "keep" in ABC approximation 
d = [5 10]; % test design

J=[5 10 20 30 40 50 60 70 80 90 100]; % number of Monte Carlo draws to test for
repeats = 100;

for j = 1:length(J)
    J(j)
	
	% MC Integration with standard pseudo random numbers
    u = zeros(repeats,1);
    for r = 1:repeats
        u(r) = utility_discrim_regen(d,J(j),prior_sims,alpha);
    end
    u_std(j) = std(u);
    u_mean(j) = mean(u);
    
	% MC Integration with randomised quasi Monte Carlo random numbers
    for r = 1:repeats
        u(r) = utility_discrim_rqmc_regen(d,J(j),prior_sims,alpha);
    end
    u_stdr(j) = std(u);
    u_meanr(j) = mean(u);
    
end

% plot results
plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')





