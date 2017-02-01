
%% Compiling of C code (requires that Matlab has access to a suitable C compiler)
mex simulate_SI.c;

%% ESTIMATION of SI model

prior_sims = 100000; % number of ABC importance samples
alpha = 0.01; % proportion of samples to "keep" in ABC approximation 
d = [1 2 3 4];  % test design

% generate from prior predictive of SI model
r = rand(prior_sims,2+50); % generate pseudo random numbers
theta_sim(:,1) = exp(norminv(r(:,1),-3.5,sqrt(0.1024)));
theta_sim(:,2) = exp(norminv(r(:,2),-4.5,sqrt(0.16)));
ysim=zeros(prior_sims,length(d));
for j = 1:prior_sims
    ysim(j,:) = simulate_SI(d,theta_sim(j,1),theta_sim(j,2),50,r(j,3:end));
end

std_sim = std(ysim);
std_sim_rep = repmat(std_sim,prior_sims,1);

J=[5 10 20 50 70 100]; % number of Monte Carlo draws to test for
repeats = 50;

% loops can be performed using parfor (parallel computing) if desired
for j = 1:length(J)
    J(j)
    u = zeros(repeats,1);
    for r = 1:repeats
        u(r) = utility_si(d,J(j),theta_sim,prior_sims,ysim,std_sim_rep,alpha);
    end
    u_std(j) = std(u);
    
    for r = 1:repeats
        u(r) = utility_si_rqmc(d,J(j),theta_sim,prior_sims,ysim,std_sim_rep,alpha);
    end
    u_stdr(j) = std(u);
    
end

% plot results
plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')



%% ESTIMATION of SI model (regeneration of prior predictive data)

J=[5 10 20 50 70 100]; % number of Monte Carlo draws to test for
repeats = 50;
prior_sims = 100000;
alpha = 0.001;
d = [1 2 3 4];

% loops can be performed using parfor (parallel computing) if desired
for j = 1:length(J)
    J(j)
    u = zeros(repeats,1);
    ur = zeros(repeats,1);

    for r = 1:repeats
        u(r) = utility_si_regen(d,J(j),prior_sims,alpha);
    end
    u_std(j) = std(u);
    
    for r = 1:repeats
        ur(r) = utility_si_regen_rqmc(d,J(j),prior_sims,alpha);
    end
    u_stdr(j) = std(ur);

end

% plot results
plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')





