% Analysis of the baseline model
%
%% Table 2: simulations baseline
clear;close all;clc;
load('simulations_baseline.mat')
disp('   Default   Spreads  Debt service   B/Y      I/Y     G/Y')
disp([mean(ED) mean(ES) mean(EBservice) mean(EB)/mean(EY) mean(EI)/mean(EY) mean(EG)/mean(EY)])

disp('sd output sd i sd c sd g')
disp([mean(SDY) mean(SDI) mean(SDC) mean(SDG)])

%% Figure 5: perturbation exercise 

clear;close all;clc;
load('simulations_baseline.mat')

T = 40;
Nsim = 1000;
% Allocate memory
irf_good.z_sim = zeros(T,Nsim);
irf_good.b_sim = zeros(T,Nsim);
irf_good.g_sim = zeros(T,Nsim);
irf_good.k_sim = zeros(T,Nsim);
irf_good.default_sim = zeros(T,Nsim);
irf_good.default_decision_sim = zeros(T,Nsim);
irf_good.q_sim = zeros(T,Nsim);
irf_good.l_sim = zeros(T,Nsim);
irf_good.y_sim = zeros(T,Nsim);

irf_bad.z_sim = zeros(T,Nsim);
irf_bad.b_sim = zeros(T,Nsim);
irf_bad.g_sim = zeros(T,Nsim);
irf_bad.k_sim = zeros(T,Nsim);
irf_bad.default_sim = zeros(T,Nsim);
irf_bad.default_decision_sim = zeros(T,Nsim);
irf_bad.q_sim = zeros(T,Nsim);
irf_bad.l_sim = zeros(T,Nsim);
irf_bad.y_sim = zeros(T,Nsim);

irf_very_bad.z_sim = zeros(T,Nsim);
irf_very_bad.b_sim = zeros(T,Nsim);
irf_very_bad.g_sim = zeros(T,Nsim);
irf_very_bad.k_sim = zeros(T,Nsim);
irf_very_bad.default_sim = zeros(T,Nsim);
irf_very_bad.default_decision_sim = zeros(T,Nsim);
irf_very_bad.q_sim = zeros(T,Nsim);
irf_very_bad.l_sim = zeros(T,Nsim);
irf_very_bad.y_sim = zeros(T,Nsim);


% steady state values
ss.y = mean(EY);
ss.z = mean(EZ);
ss.k = mean(EK);
ss.b = mean(EB);
ss.g = mean(EG);
ss.i = mean(EI);
ss.q = mean(EQ);

% Simulations good
zero_b_index = 1;
b_init = interp1(grid.b,1:N.b, mean(EB), 'nearest', 'extrap');
k_init = interp1(grid.k,1:N.k, mean(EK), 'nearest', 'extrap');
z_init = interp1(grid.z, 1:N.z, mean(grid.z)+3*sig, 'nearest');
X0 = zeros(1,N.z);
X0(z_init) = 1;


for j=1:Nsim
    disp('Simulation no')
    display(num2str(j))
	% Markov Chain simulation
	mc = dtmc(P);
	z_sim_indices = simulate(mc, T+1,'X0',X0);

	% Initial state
    in_default = 0; 
	ib = b_init;
	ik = k_init;
    
    for t=1:T
		iz = z_sim_indices(t);
  
        irf_good.b_sim(t,j) = grid.b(ib);
        irf_good.k_sim(t,j) = grid.k(ik);
  		
        % Decision tree
        if in_default==0
            if default_states(iz,ib,ik)>eps
                in_default=1;
                default_decision=1;
                ibnext = zero_b_index;
                iknext = p.ik_d(iz,ik);
                irf_good.g_sim(t,j) = p.g_d(iz,ik);
            else
                default_decision=0;
                ibnext = p.ib(iz,ib,ik);
                iknext = p.ik_c(iz,ib,ik);
                irf_good.g_sim(t,j) = p.g_c(iz,ib,ik);
            end
        else
            default_decision=0;
            ibnext = zero_b_index;
            iknext = p.ik_d(iz,ik);
            irf_good.g_sim(t,j) = p.g_d(iz,ik);
            if rand<theta
                in_default=0;
            end
        end
        
        irf_good.default_sim(t,j) = in_default;
        irf_good.default_decision_sim(t,j) = default_decision;
        if in_default == 0
            irf_good.z_sim(t,j) = grid.z(iz);
        else
            irf_good.z_sim(t,j) = min(grid.z(iz), phi*mean(grid.z));
        end
        irf_good.q_sim(t,j) = q(iz,ibnext,iknext);

        ib = ibnext;
		ik = iknext;
    end
   
    % Static variables 
    irf_good.l_sim(:,j) = (((1-alpha).*zeta*irf_good.z_sim(:,j).*irf_good.k_sim(:,j).^alpha)./((1+tau))).^(1/(psi+alpha));
    irf_good.y_sim(:,j) = zeta*irf_good.z_sim(:,j).*(irf_good.k_sim(:,j).^alpha).*irf_good.l_sim(:,j).^(1-alpha);
end

irf_good.decision_d = mean(irf_good.default_decision_sim,2);
irf_good.y = (mean(irf_good.y_sim,2)-ss.y)/ss.y;
irf_good.k = (mean(irf_good.k_sim,2)-ss.k)/ss.k;
irf_good.z = (mean(irf_good.z_sim,2)-ss.z)/ss.z;
irf_good.g = (mean(irf_good.g_sim,2)-ss.g)/ss.g;
irf_good.b = (mean(irf_good.b_sim,2)-ss.b)/ss.b;
irf_good.q = (mean(irf_good.q_sim,2)-ss.q)/ss.q;
irf_good.i = zeros(T,1);
irf_good.i(2:end) = irf_good.k(2:end) - (1-delta)*irf_good.k(1:end-1);

% Simulations bad
z_init = interp1(grid.z, 1:N.z, mean(grid.z)-3*sig, 'nearest');
X0 = zeros(1,N.z);
X0(z_init) = 1;

for j=1:Nsim
    disp('Simulation no')
    display(num2str(j))
	% Markov Chain simulation
	mc = dtmc(P);
	z_sim_indices = simulate(mc, T+1,'X0',X0);

	% Initial state
    in_default = 0; 
	ib = b_init;
	ik = k_init;
    
    for t=1:T

        iz = z_sim_indices(t);

        irf_bad.b_sim(t,j) = grid.b(ib);
        irf_bad.k_sim(t,j) = grid.k(ik);
  		
        % Decision tree
        if in_default==0
            if default_states(iz,ib,ik)>eps
                in_default=1;
                default_decision=1;
                ibnext = zero_b_index;
                iknext = p.ik_d(iz,ik);
                irf_bad.g_sim(t,j) = p.g_d(iz,ik);
            else
                default_decision=0;
                ibnext = p.ib(iz,ib,ik);
                iknext = p.ik_c(iz,ib,ik);
                irf_bad.g_sim(t,j) = p.g_c(iz,ib,ik);
            end
        else
            default_decision=0;
            ibnext = zero_b_index;
            iknext = p.ik_d(iz,ik);
            irf_bad.g_sim(t,j) = p.g_d(iz,ik);
            if rand<theta
                in_default=0;
            end
        end
        
        irf_bad.default_sim(t,j) = in_default;
        irf_bad.default_decision_sim(t,j) = default_decision;
        if in_default == 0
            irf_bad.z_sim(t,j) = grid.z(iz);
        else
            irf_bad.z_sim(t,j) = min(grid.z(iz), phi*mean(grid.z));
        end
        irf_bad.q_sim(t,j) = q(iz,ibnext,iknext);
        
		ib = ibnext;
		ik = iknext;
    end
    % Static variables
    irf_bad.l_sim(:,j) = (((1-alpha).*zeta*irf_bad.z_sim(:,j).*irf_bad.k_sim(:,j).^alpha)./((1+tau))).^(1/(psi+alpha));
    irf_bad.y_sim(:,j) = zeta*irf_bad.z_sim(:,j).*(irf_bad.k_sim(:,j).^alpha).*irf_bad.l_sim(:,j).^(1-alpha);
end

irf_bad.decision_d = mean(irf_bad.default_decision_sim,2);
irf_bad.y = (mean(irf_bad.y_sim,2)-ss.y)/ss.y;
irf_bad.k = (mean(irf_bad.k_sim,2)-ss.k)/ss.k;
irf_bad.z = (mean(irf_bad.z_sim,2)-ss.z)/ss.z;
irf_bad.g = (mean(irf_bad.g_sim,2)-ss.g)/ss.g;
irf_bad.b = (mean(irf_bad.b_sim,2)-ss.b)/ss.b;
irf_bad.q = (mean(irf_bad.q_sim,2)-ss.q)/ss.q;
irf_bad.i = zeros(T,1);
irf_bad.i(2:end) = irf_bad.k(2:end) - (1-delta)*irf_bad.k(1:end-1);

%%%%
% Simulations very bad
z_init = interp1(grid.z, 1:N.z, mean(grid.z)-5*sig, 'nearest');
X0 = zeros(1,N.z);
X0(z_init) = 1;


for j=1:Nsim
    disp('Simulation no')
    display(num2str(j))
	% Markov Chain simulation
	mc = dtmc(P);
	z_sim_indices = simulate(mc, T+1,'X0',X0);

	% Initial state
    in_default = 0; 
	ib = b_init;
	ik = k_init;
    
    for t=1:T

        iz = z_sim_indices(t);

        irf_very_bad.b_sim(t,j) = grid.b(ib);
        irf_very_bad.k_sim(t,j) = grid.k(ik);
  		
        % Decision tree
        if in_default==0
            if default_states(iz,ib,ik)>eps
                in_default=1;
                default_decision=1;
                ibnext = zero_b_index;
                iknext = p.ik_d(iz,ik);
                irf_very_bad.g_sim(t,j) = p.g_d(iz,ik);
            else
                default_decision=0;
                ibnext = p.ib(iz,ib,ik);
                iknext = p.ik_c(iz,ib,ik);
                irf_very_bad.g_sim(t,j) = p.g_c(iz,ib,ik);
            end
        else
            default_decision=0;
            ibnext = zero_b_index;
            iknext = p.ik_d(iz,ik);
            irf_very_bad.g_sim(t,j) = p.g_d(iz,ik);
            if rand<theta
                in_default=0;
            end
        end
        
        irf_very_bad.default_sim(t,j) = in_default;
        irf_very_bad.default_decision_sim(t,j) = default_decision;
        if in_default == 0
            irf_very_bad.z_sim(t,j) = grid.z(iz);
        else
            irf_very_bad.z_sim(t,j) = min(grid.z(iz), phi*mean(grid.z));
        end
        irf_very_bad.q_sim(t,j) = q(iz,ibnext,iknext);

        ib = ibnext;
		ik = iknext;
    end
   
    irf_very_bad.l_sim(:,j) = (((1-alpha).*zeta*irf_very_bad.z_sim(:,j).*irf_very_bad.k_sim(:,j).^alpha)./((1+tau))).^(1/(psi+alpha));
    irf_very_bad.y_sim(:,j) = zeta*irf_very_bad.z_sim(:,j).*(irf_very_bad.k_sim(:,j).^alpha).*irf_very_bad.l_sim(:,j).^(1-alpha);
end

irf_very_bad.decision_d = mean(irf_very_bad.default_decision_sim,2);
irf_very_bad.y = (mean(irf_very_bad.y_sim,2)-ss.y)/ss.y;
irf_very_bad.k = (mean(irf_very_bad.k_sim,2)-ss.k)/ss.k;
irf_very_bad.z = (mean(irf_very_bad.z_sim,2)-ss.z)/ss.z;
irf_very_bad.g = (mean(irf_very_bad.g_sim,2)-ss.g)/ss.g;
irf_very_bad.b = (mean(irf_very_bad.b_sim,2)-ss.b)/ss.b;
irf_very_bad.q = (mean(irf_very_bad.q_sim,2)-ss.q)/ss.q;
irf_very_bad.i = zeros(T,1);
irf_very_bad.i(2:end) = irf_very_bad.k(2:end) - (1-delta)*irf_very_bad.k(1:end-1);

newcolors = [0.00 0.20 1.00
             1.00 0.54 0.00
             0.83 0.14 0.14];
         
colororder(newcolors)

figure(1)
subplot(4,2,1), plot([irf_good.y irf_bad.y irf_very_bad.y], 'Linewidth', 1.2)
title('Output','interpreter', 'latex', 'fontsize', 12)
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,2), plot([irf_good.z irf_bad.z irf_very_bad.z], 'Linewidth', 1.2)
title('Productivity','interpreter', 'latex', 'fontsize', 12)
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,3), plot([irf_good.g irf_bad.g irf_very_bad.g ], 'Linewidth', 1.2)
title('Public expenditure','interpreter', 'latex', 'fontsize', 12)
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,4), plot([irf_good.k(2:end) irf_bad.k(2:end) irf_very_bad.k(2:end)], 'Linewidth', 1.2)
title('Public capital','interpreter', 'latex', 'fontsize', 12), xticks([1 5 10 15 20 25 30 35 40])
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,5),plot([irf_good.i(2:end) irf_bad.i(2:end) irf_very_bad.i(2:end)], 'Linewidth', 1.2)
title('Public investment','interpreter', 'latex', 'fontsize', 12), xticks([1 5 10 15 20 25 30 35 40])
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,6), plot([irf_good.b(2:end) irf_bad.b(2:end) irf_very_bad.b(2:end)], 'Linewidth', 1.2)
title('Bond stock','interpreter', 'latex', 'fontsize', 12), xticks([1 5 10 15 20 25 30 35 40])
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,7), plot([irf_good.q irf_bad.q ], 'Linewidth', 1.2)
title('Bond prices','interpreter', 'latex', 'fontsize', 12)
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Dev.', 'interpreter', 'latex')
subplot(4,2,8), plot([irf_good.decision_d irf_bad.decision_d irf_very_bad.decision_d], 'Linewidth', 1.2)
title('Default','interpreter', 'latex', 'fontsize', 12)
xlabel('Quarters', 'interpreter', 'latex'), ylabel('\% Freq.', 'interpreter', 'latex')
legend1 = legend('show', '$+3\sigma_{\varepsilon}$', '$-3\sigma_{\varepsilon}$','$-5\sigma_{\varepsilon}$' );
set(legend1,...
    'Position',[0.358 0.025 0.354 0.043],...
    'NumColumns',3,...
    'Interpreter','latex',...
    'FontAngle','italic',...
    'FontSize',11,...
    'FontName','Bell MT');
title(legend1,'Initial shocks');

%% Figure 6: default episodes 
clear;close all;clc;
load('simulations_baseline.mat')
episodes = 16;
count = 0;
for j = 1:Nsim
    for t = 1:T-burn-episodes
        if default_decision_sim(t,j)==1 & t > episodes 
            count = count +1;
            def.y(count,1:2*episodes+1) = y_sim(t-episodes:t+episodes,j);
            def.b(count,1:2*episodes+1) = b_sim(t-episodes:t+episodes,j);
            def.k(count,1:2*episodes+1) = k_sim(t-episodes:t+episodes,j);
            def.g(count,1:2*episodes+1) = g_sim(t-episodes:t+episodes,j);
            def.q(count,1:episodes) = q_sim(t-episodes:t-1,j);
            def.i(count,1:2*episodes+1) = i_sim(t-episodes:t+episodes,j);
        end
    end
end
def.y = mean(def.y);
def.b = mean(def.b);
def.k = mean(def.k);
def.g = mean(def.g);
def.q = mean(def.q);
def.i = mean(def.i);

def.gy = def.g./def.y;
def.by = def.b./def.y;
def.ky = def.k./def.y;


figure(1)
subplot(3,2,1), plot(def.y, 'Linewidth', 1.2),  xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes);
title('$y_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)
subplot(3,2,2), plot([def.g' def.gy'], 'Linewidth', 1.2),  xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes), legend('show', 'g', 'g/y','interpreter', 'latex', 'fontsize', 11,'location', 'northwest');
title('$g_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)
subplot(3,2,3), plot([def.b' def.by'], 'Linewidth', 1.2),  xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes),legend('show', 'b', 'b/y','interpreter', 'latex', 'fontsize', 11, 'location', 'southwest');
title('$b_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)
subplot(3,2,4), plot(def.q, 'Linewidth', 1.2),  xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes);
title('$q_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)
subplot(3,2,5), plot([def.k' def.ky'], 'Linewidth', 1.2),  xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes),legend('show', 'k', 'k/y','interpreter', 'latex', 'fontsize', 11,'location', 'northwest');
title('$k_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)
subplot(3,2,6), plot(def.i(2:end), 'Linewidth', 1.2), xline(episodes+1, 'Linewidth', 1.2);
xticks(1:4:2*episodes+1), xticklabels(-episodes:4:episodes);
title('$i_t$','interpreter', 'latex', 'fontsize', 12), xlabel('Quarters','interpreter', 'latex', 'fontsize', 12)







