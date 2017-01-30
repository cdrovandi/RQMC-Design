
% parpool(4);    % parallel computing with 4 cores
warning off all;
% pctRunOnAll warning('off','all');   % use this command if using parfor so that the warnings are turned off on all cores
 
d = 1:1:15; % test design

J=[5 10 20 30 40 50 60 70 80 90 100]; % number of Monte Carlo draws to test for
repeats = 50;

u_std = zeros(length(J),1);
u_stdr = zeros(length(J),1);

% loops can be performed using parfor (parallel computing) if desired
for j = 1:length(J)
    J(j)
    u = zeros(repeats,1);
    ur = zeros(repeats,1);
    for r = 1:repeats
        u(r) = utility_pk(d,J(j));
    end
    u_std(j) = std(u);
    
    for r = 1:repeats
        ur(r) = utility_pk_rqmc(d,J(j));
    end
    u_stdr(j) = std(ur);
    
end

plot(J,u_std,'o--','LineWidth',2)
hold on;
plot(J,u_stdr,'x--','LineWidth',2)
xlabel('J','FontSize',16);
ylabel('sd($$U_J(d)$$)','Interpreter','Latex','FontSize',16);
set(gca,'XTick',J);
set(gca, 'YTick', 0:0.1:50);
legend('MC','RQMC')



