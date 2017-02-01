
%% LOGISTIC REGRESSION EXAMPLE fixed prior simulations
load('Xrand.mat'); % test design

% prior specifications
mu_prior = [0,7,8,-3,0.5];
std_prior = sqrt([3,3,3,3,3]);

sims = 100000; % number of samples for importance sampling

% draw from prior
beta_prior = normrnd(repmat(mu_prior,sims,1), repmat(std_prior,sims,1));

J=[10 20 50 100]; % number of Monte Carlo draws to test for
repeats = 100;

u_std = zeros(length(J),1);
u_stdr = zeros(length(J),1);

% loops can be performed using parfor (parallel computing) if desired
for i = 1:length(J)
    J(i)
    u = zeros(repeats,1);
    ur = zeros(repeats,1);
    for j = 1:repeats
        u(j) = utility_logreg(X,J(i),beta_prior,mu_prior,std_prior);
    end
    u_std(i) = std(u);
    
    for j = 1:repeats
        ur(j) = utility_logreg_rqmc(X,J(i),beta_prior,mu_prior,std_prior);
    end
    u_stdr(i) = std(ur);
    
end

plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')


%% LOGISTIC REGRESSION EXAMPLE regenerated prior simulations

load('Xrand.mat'); % test design

% prior specifications
mu_prior = [0,7,8,-3,0.5];
std_prior = sqrt([3,3,3,3,3]);

sims = 100000; % number of samples for importance sampling

J=[10 20 50 100]; % number of Monte Carlo draws to test for
repeats = 100;
u = zeros(repeats,1);
ur = zeros(repeats,1);

u_std = zeros(length(J),1);
u_stdr = zeros(length(J),1);

% loops can be performed using parfor (parallel computing) if desired
for i = 1:length(J)
    J(i)
    for j = 1:repeats
        u(j) = utility_logreg_est(X,J(i),sims,mu_prior,std_prior);
    end
    u_std(i) = std(u);
    
    for j = 1:repeats
        ur(j) = utility_logreg_est_rqmc(X,J(i),sims,mu_prior,std_prior);
    end
    u_stdr(i) = std(ur);
    
end

plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
legend('MC','RQMC')


