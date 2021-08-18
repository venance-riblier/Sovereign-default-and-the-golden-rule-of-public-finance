clear;close all;clc;
load('policies.mat')

Nsim = 1000;
T = 4500;
burn= 500;
    
% Allocate memory
z_sim = zeros(T,Nsim);
b_sim = zeros(T,Nsim);
g_sim = zeros(T,Nsim);
k_sim = zeros(T,Nsim);
default_sim = zeros(T,Nsim);
default_decision_sim = zeros(T,Nsim);
q_sim = zeros(T,Nsim);
b_service = zeros(T,Nsim);
b_y = zeros(T,Nsim);
l_sim = zeros(T,Nsim);
y_sim = zeros(T,Nsim);
c_sim = zeros(T,Nsim);
tax_sim = zeros(T,Nsim);
i_sim = zeros(T,Nsim);
spread = zeros(T,Nsim);
u_sim = zeros(T,Nsim);



% Set initial conditions
zero_b_index = 1;
z_init = (N.z+1)/2;
X0 = zeros(1,N.z);
X0(z_init) = randi(N.z,1);


for j=1:Nsim
    disp('Simulation no')
    display(num2str(j))
	% Markov Chain simulation
	mc = dtmc(P);
	z_sim_indices = simulate(mc, T+1,'X0',X0);
   
	% Initial state
    in_default = 0; 
	ib = randi(N.b,1); 
	ik = randi(N.k,1);
    
    for t=1:T
        % Random state
		iz = z_sim_indices(t);
        % Store variable
        b_sim(t,j) = grid.b(ib);
        k_sim(t,j) = grid.k(ik);
  		
        % Decision tree
        if in_default==0
            if default_states(iz,ib,ik)>eps
                in_default=1;
                default_decision=1;
                ibnext = zero_b_index;
                iknext = p.ik_d(iz,ik);
                g_sim(t,j) = p.g_d(iz,ik);
            else
                default_decision=0;
                ibnext = p.ib(iz,ib,ik);
                iknext = p.ik_c(iz,ib,ik);
                g_sim(t,j) = p.g_c(iz,ib,ik);
            end
        else
            default_decision=0;
            ibnext = zero_b_index;
            iknext = p.ik_d(iz,ik);
            g_sim(t,j) = p.g_d(iz,ik);
            if rand<theta
                in_default=0;
            end
        end
        % Store variables
        default_sim(t,j) = in_default;
        default_decision_sim(t,j) = default_decision;
        if in_default == 0
            z_sim(t,j) = grid.z(iz);
        else
            z_sim(t,j) = min(grid.z(iz), phi*mean(grid.z));
        end
        q_sim(t,j) = q(iz,ibnext,iknext);
        % Actualize State Space
		ib = ibnext;
		ik = iknext;
    end
    
    l_sim(:,j) = (((1-alpha).*zeta*z_sim(:,j).*k_sim(:,j).^alpha)./((1+tau))).^(1/(psi+alpha));
    y_sim(:,j) = zeta*z_sim(:,j).*(k_sim(:,j).^alpha).*l_sim(:,j).^(1-alpha);
    b_y(:,j) = b_sim(:,j)./y_sim(:,j);
    b_service(:,j) = (lambda + coupon*(1-lambda))*b_sim(:,j);
    c_sim(:,j) = y_sim(:,j)./(1+tau);
    tax_sim(:,j) = tau*c_sim(:,j);
    i_sim(2:end,j) = k_sim(2:end,j) - (1-delta).*k_sim(1:end-1,j);
    rr = (lambda +(1-lambda)*coupon)./q_sim(:,j) - lambda;
    spread(:,j) = (1+rr).^4-(1+r).^4;
    u_sim(:,j) = (1-omega).*(c_sim(:,j) - l_sim(:,j).^(1+psi)./(1+psi)).^(1-sigma)./(1-sigma) + ...
        omega.*g_sim(:,j).^(1-sigma)./(1-sigma);
  
end

    
z_sim(1:burn,:) = [];
k_sim(1:burn,:) = [];
b_sim(1:burn,:) = [];
g_sim(1:burn,:) = [];
default_sim(1:burn,:) = [];
default_decision_sim(1:burn,:) = [];
q_sim(1:burn,:) = [];
b_service(1:burn,:) = [];
b_y(1:burn,:) = [];
l_sim(1:burn,:) = [];
y_sim(1:burn,:) = [];
c_sim(1:burn,:) = [];
tax_sim(1:burn,:) = [];
i_sim(1:burn,:) = [];
spread(1:burn,:) = [];
u_sim(1:burn,:) = [];

% exclude default periods
y_nd = y_sim;
for j = 1:Nsim
    x = find(default_sim(:,j)==1);
    b_service(x,j) = NaN;
    b_y(x,j) = NaN;
    y_nd(x,j) = NaN;   
end
b_service(isnan(b_service)) = [];
b_y(isnan(b_y)) = [];
y_nd(isnan(y_nd)) = [];

%% HP filter
[ZHPT, ZHP] = hpfilter(z_sim, 1600);
[YHPT, YHP] = hpfilter(y_sim, 1600);
[YNDHPT, YNDHP] = hpfilter(y_nd, 1600);
[BHPT, BHP] = hpfilter(b_sim, 1600);
[BserviceHPT, BserviceHP] = hpfilter(b_service, 1600);
[B_YHPT, B_YHP] = hpfilter(b_y, 1600);
[LHPT, LHP] = hpfilter(l_sim, 1600);
[KHPT, KHP] = hpfilter(k_sim, 1600);
[IHPT, IHP] = hpfilter(i_sim, 1600);
[GHPT, GHP] = hpfilter(g_sim, 1600);
[CHPT, CHP] = hpfilter(c_sim, 1600);
[QHPT, QHP] = hpfilter(q_sim, 1600);
[SHPT, SHP] = hpfilter(spread, 1600);
[TAXHPT, TAXHP] = hpfilter(tax_sim, 1600);
[UHPT, UHP] = hpfilter(u_sim, 1600);


% Averages values of each simulations
EY = mean(YHPT);
EZ = mean(ZHPT);
EB = mean(BHPT);
EBY = mean(B_YHPT);
EBservice = mean(BserviceHPT);
EL = mean(LHPT);
EK = mean(KHPT);
EI = mean(IHPT);
EG = mean(GHPT);
EC = mean(CHPT);
ED = sum(default_decision_sim)./((T-burn)/400);
ES = mean(SHPT);
ETAX = mean(TAXHPT);
EQ = mean(QHPT);
EU = mean(UHPT);
% Standard deviation of each simulations 
SDY = std(YNDHP);
SDC = std(CHP);
SDG = std(GHP);
SDS = std(SHP);
SDI = std(IHP);
SDB = std(BHP);
SDU = std(UHP);

%%%% Correlations %%%%
for j=1:Nsim
    
    aux=corrcoef(BHP(:,j),YHP(:,j));
    cor.BY(j) = aux(2,1);
    aux=corrcoef(GHP(:,j),YHP(:,j));
    cor.GY(j) = aux(2,1);
    aux=corrcoef(IHP(:,j),YHP(:,j));
    cor.IY(j) = aux(2,1);
    
    aux=corrcoef(CHP(:,j),YHP(:,j));
    cor.CY(j) = aux(2,1);
    
end

save('simulations_baseline.mat')

