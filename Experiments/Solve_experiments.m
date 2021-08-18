%% Golden rule experiment, ro=0
clear;close all;clc;

    % Calibration 
    rollover = 0.01;
    r = 0.01; 
    beta = 0.86; 
    theta = 0.1; 
    phi = 0.86; 

    sigma = 2; 
    psi = 1/2.22; 
    omega = 0.3; 

    alpha = 0.36; 

    rho = 0.95; 
    sig = 0.01; 

    tau = 0.18; 
    delta = 0.025; 
    lambda = 0.05; 
    coupon=0.03; 
    zeta = 1.22; 

    % State space
    % Grid points
    N.z = 9;
    N.k = 20;
    N.b = 20;

    % Debt Grid
    bmin = 0;
    bmax = 1.1;
    grid.b = linspace(bmin, bmax, N.b); 
    index0 = 1;

    % Capital Grid
    kmin = 0.15;
    kmax = 3; 
    grid.k = linspace(kmin, kmax, N.k);
    
    % TFP
    w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
    [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
    grid.z = exp(grid.z)';
    z_mean_index = (N.z+1)/2;

    
    % Initial grids
    Vd = ones(N.z,N.k);
    Vc = ones(N.z, N.b, N.k);
    V = Vc;
    V1 = Vc;
    default_states = zeros(N.z, N.b, N.k); 
    q =  ones(N.z, N.b, N.k)*(lambda+(1-lambda)*coupon)/(r+lambda);

    
    % VFI loop
    dist = 'Start';
    iter = 0;
    tol = 1e-6;
    while dist>tol        
        iter = iter+1;
        display(num2str([iter dist]))
        
        % Value of default 
        for ik = 1:N.k
            for iz = 1:N.z
                k_d = grid.k(ik);
                z_d = min(grid.z(iz), phi*mean(grid.z));
                l_d = (((1-alpha)*z_d*zeta*k_d^alpha)/((1+tau)))^(1/(psi+alpha));
                y_d = zeta*z_d*k_d^alpha*l_d^(1-alpha);
                C_d =  y_d/(1+tau);
                i_d =  grid.k - (1-delta)*k_d;
                G_d = tau*C_d - i_d;
                neg = G_d<=0;
                constraint_d = G_d > tau*C_d;
                u_d = (C_d -l_d^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                V_d_temp = (1-omega)*u_d +  omega*G_d.^(1-sigma)/(1-sigma) -neg*10^16 - constraint_d*10^16;

                % Compute expectations
                for iznext = 1:N.z
                    aux = reshape(V(iznext,index0,:), 1, N.k);
                    V_d_temp = V_d_temp + beta*P(iz,iznext).*(theta*aux + (1-theta)*Vd(iznext,:)); 
                end
                [Vd(iz,ik), ik_d_opt] = max(V_d_temp);
                p.ik_d(iz,ik) = ik_d_opt;
                p.k_d(iz,ik) = grid.k(ik_d_opt);
                p.g_d(iz,ik) = G_d(ik_d_opt);
            end
        end

        % Value of staying in debt contract 
        for ib = 1:N.b
            for ik = 1:N.k
                for iz = 1:N.z
                    b = grid.b(ib);
                    qq = q(iz,:,:);
                    qq = reshape(qq, [N.b, N.k]);
                    l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                    C = y_c/(1+tau);
                    nfa = qq.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                    i = grid.k - (1-delta)*grid.k(ik);
                    G = tau*C - i + nfa;
                    neg2 = G<=0  ;
                    constraint_c = G + (1-rollover)*(lambda+(1-lambda)*coupon)*grid.b(ib) > tau*C;
                    uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16 - constraint_c*10^16; 
             
                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,:,:), N.b, N.k);
                        V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                    end            
                    
                    % Compute policy function
                    [Vc(iz,ib,ik), pol] = maxmat(V_c_temp);
                    p.ib(iz,ib,ik) = pol(1);
                    p.ik_c(iz,ib,ik) = pol(2);
                    p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                    p.k_c(iz,ib,ik) = grid.k(p.ik_c(iz,ib,ik));
                    p.g_c(iz,ib,ik) = G(pol(1), pol(2));
                                
                    V1(iz,ib,ik) = max(Vc(iz,ib,ik),Vd(iz,ik));
                    default_states(iz,ib,ik) = 1*(Vd(iz,ik)>Vc(iz,ib,ik));
                    
                    if default_states(iz,ib,ik)==1
                        p.ib(iz,ib,ik) = index0;
                        p.b(iz,ib,ik) = grid.b(index0);
                        p.k(iz,ib,ik) = p.k_d(iz,ik);
                        p.ik(iz,ib,ik) = p.ik_d(iz,ik);
                        p.g(iz,ib,ik) = p.g_d(iz,ik);
                    else
                        p.k(iz,ib,ik) = p.k_c(iz,ib,ik);
                        p.ik(iz,ib,ik) = p.ik_c(iz,ib,ik);
                        p.g(iz,ib,ik) = p.g_c(iz,ib,ik);
                    end
                end
            end
        end
        dist = max(max(max(abs(V1-V))));
        V = V1;
    
        % Update q
        q1 = ones(N.z,N.b,N.k);
        for ik = 1:N.k
            repay = 1-default_states(:,:,ik);
            q1(:,:,ik) = P*repay.*(lambda + (1-lambda)*(coupon+q(:,:,ik)));
        end
        q = q1./(1+r);
        
    end
    
    save('Policies\GR_ro_0.mat')
    
%% GR experiment, ro=0.5
clear;close all;clc;

    % Calibration 
    rollover = 0.5;
    r = 0.01; 
    beta = 0.86; 
    theta = 0.1; 
    phi = 0.86; 

    sigma = 2; 
    psi = 1/2.22; 
    omega = 0.3; 

    alpha = 0.36; 

    rho = 0.95; 
    sig = 0.01; 

    tau = 0.18; 
    delta = 0.025; 
    lambda = 0.05; 
    coupon=0.03; 
    zeta = 1.22; 

    % State space
    % Grid points
    N.z = 9;
    N.k = 20;
    N.b = 30;

    % Debt Grid
    bmin = 0;
    bmax = 0.9;
    grid.b = linspace(bmin, bmax, N.b); 
    index0 = 1;

    % Capital Grid
    kmin = 0.15;
    kmax = 3; 
    grid.k = linspace(kmin, kmax, N.k);
    
    % TFP
    w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
    [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
    grid.z = exp(grid.z)';
    z_mean_index = (N.z+1)/2;

    
    % Initial grids
    Vd = ones(N.z,N.k);
    Vc = ones(N.z, N.b, N.k);
    V = Vc;
    V1 = Vc;
    default_states = zeros(N.z, N.b, N.k); 
    q =  ones(N.z, N.b, N.k)*(lambda+(1-lambda)*coupon)/(r+lambda);

    
    % VFI loop
    dist = 'Start';
    iter = 0;
    tol = 1e-6;
    while dist>tol        
        iter = iter+1;
        display(num2str([iter dist]))
        
        % Value of default 
        for ik = 1:N.k
            for iz = 1:N.z
                k_d = grid.k(ik);
                z_d = min(grid.z(iz), phi*mean(grid.z));
                l_d = (((1-alpha)*z_d*zeta*k_d^alpha)/((1+tau)))^(1/(psi+alpha));
                y_d = zeta*z_d*k_d^alpha*l_d^(1-alpha);
                C_d =  y_d/(1+tau);
                i_d =  grid.k - (1-delta)*k_d;
                G_d = tau*C_d - i_d;
                neg = G_d<=0;
                constraint_d = G_d > tau*C_d;
                u_d = (C_d -l_d^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                V_d_temp = (1-omega)*u_d +  omega*G_d.^(1-sigma)/(1-sigma) -neg*10^16 - constraint_d*10^16;

                % Compute expectations
                for iznext = 1:N.z
                    aux = reshape(V(iznext,index0,:), 1, N.k);
                    V_d_temp = V_d_temp + beta*P(iz,iznext).*(theta*aux + (1-theta)*Vd(iznext,:)); 
                end
                [Vd(iz,ik), ik_d_opt] = max(V_d_temp);
                p.ik_d(iz,ik) = ik_d_opt;
                p.k_d(iz,ik) = grid.k(ik_d_opt);
                p.g_d(iz,ik) = G_d(ik_d_opt);
            end
        end
 
 
        
        % Value of staying in debt contract 
        for ib = 1:N.b
            for ik = 1:N.k
                for iz = 1:N.z
                    b = grid.b(ib);
                    qq = q(iz,:,:);
                    qq = reshape(qq, [N.b, N.k]);
                    l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                    C = y_c/(1+tau);
                    nfa = qq.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                    i = grid.k - (1-delta)*grid.k(ik);
                    G = tau*C - i + nfa;
                    neg2 = G<=0  ;
                    constraint_c = G + (1-rollover)*(lambda+(1-lambda)*coupon)*grid.b(ib) > tau*C;
                    uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16 - constraint_c*10^16; 
             
                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,:,:), N.b, N.k);
                        V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                    end            
                    
                    % Compute policy function
                    [Vc(iz,ib,ik), pol] = maxmat(V_c_temp);
                    p.ib(iz,ib,ik) = pol(1);
                    p.ik_c(iz,ib,ik) = pol(2);
                    p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                    p.k_c(iz,ib,ik) = grid.k(p.ik_c(iz,ib,ik));
                    p.g_c(iz,ib,ik) = G(pol(1), pol(2));
                    V1(iz,ib,ik) = max(Vc(iz,ib,ik),Vd(iz,ik));
                    default_states(iz,ib,ik) = 1*(Vd(iz,ik)>Vc(iz,ib,ik));
                    
                    if default_states(iz,ib,ik)==1
                        p.ib(iz,ib,ik) = index0;
                        p.b(iz,ib,ik) = grid.b(index0);
                        p.k(iz,ib,ik) = p.k_d(iz,ik);
                        p.ik(iz,ib,ik) = p.ik_d(iz,ik);
                        p.g(iz,ib,ik) = p.g_d(iz,ik);
                    else
                        p.k(iz,ib,ik) = p.k_c(iz,ib,ik);
                        p.ik(iz,ib,ik) = p.ik_c(iz,ib,ik);
                        p.g(iz,ib,ik) = p.g_c(iz,ib,ik);
                    end
                end
            end
        end
        dist = max(max(max(abs(V1-V))));
        V = V1;
    
        % Update q
        q1 = ones(N.z,N.b,N.k);
        for ik = 1:N.k
            repay = 1-default_states(:,:,ik);
            q1(:,:,ik) = P*repay.*(lambda + (1-lambda)*(coupon+q(:,:,ik)));
        end
        q = q1./(1+r);
    end
    
    save('Policies\GR_ro_0.5.mat')
    
%% GR experiment, ro=1
clear;close all;clc;

    % Calibration 
    rollover = 1;
    r = 0.01; 
    beta = 0.86; 
    theta = 0.1; 
    phi = 0.86; 

    sigma = 2; 
    psi = 1/2.22; 
    omega = 0.3; 

    alpha = 0.36; 

    rho = 0.95; 
    sig = 0.01; 

    tau = 0.18; 
    delta = 0.025; 
    lambda = 0.05; 
    coupon=0.03; 
    zeta = 1.22; 

    % State space
    % Grid points
    N.z = 9;
    N.k = 20;
    N.b = 20;

    % Debt Grid
    bmin = 0;
    bmax = 1.1;
    grid.b = linspace(bmin, bmax, N.b); 
    index0 = 1;

    % Capital Grid
    kmin = 0.15;
    kmax = 3; 
    grid.k = linspace(kmin, kmax, N.k);
    
    % TFP
    w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
    [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
    grid.z = exp(grid.z)';
    z_mean_index = (N.z+1)/2;

    
    % Initial grids
    Vd = ones(N.z,N.k);
    Vc = ones(N.z, N.b, N.k);
    V = Vc;
    V1 = Vc;
    default_states = zeros(N.z, N.b, N.k); 
    q =  ones(N.z, N.b, N.k)*(lambda+(1-lambda)*coupon)/(r+lambda);

    
    % VFI loop
    dist = 'Start';
    iter = 0;
    tol = 1e-6;
    while dist>tol        
        iter = iter+1;
        display(num2str([iter dist]))
        
        % Value of default 
        for ik = 1:N.k
            for iz = 1:N.z
                k_d = grid.k(ik);
                z_d = min(grid.z(iz), phi*mean(grid.z));
                l_d = (((1-alpha)*z_d*zeta*k_d^alpha)/((1+tau)))^(1/(psi+alpha));
                y_d = zeta*z_d*k_d^alpha*l_d^(1-alpha);
                C_d =  y_d/(1+tau);
                i_d =  grid.k - (1-delta)*k_d;
                G_d = tau*C_d - i_d;
                neg = G_d<=0;
                constraint_d = G_d > tau*C_d;
                u_d = (C_d -l_d^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                V_d_temp = (1-omega)*u_d +  omega*G_d.^(1-sigma)/(1-sigma) -neg*10^16 - constraint_d*10^16;

                % Compute expectations
                for iznext = 1:N.z
                    aux = reshape(V(iznext,index0,:), 1, N.k);
                    V_d_temp = V_d_temp + beta*P(iz,iznext).*(theta*aux + (1-theta)*Vd(iznext,:)); 
                end
                [Vd(iz,ik), ik_d_opt] = max(V_d_temp);
                p.ik_d(iz,ik) = ik_d_opt;
                p.k_d(iz,ik) = grid.k(ik_d_opt);
                p.g_d(iz,ik) = G_d(ik_d_opt);
            end
        end
 
        % Value of staying in debt contract 
        for ib = 1:N.b
            for ik = 1:N.k
                for iz = 1:N.z
                    b = grid.b(ib);
                    qq = q(iz,:,:);
                    qq = reshape(qq, [N.b, N.k]);
                    l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                    C = y_c/(1+tau);
                    nfa = qq.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                    i = grid.k - (1-delta)*grid.k(ik);
                    G = tau*C - i + nfa;
                    neg2 = G<=0  ;
                    constraint_c = G + (1-rollover)*(lambda+(1-lambda)*coupon)*grid.b(ib) > tau*C;
                    uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16 - constraint_c*10^16; 
             
                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,:,:), N.b, N.k);
                        V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                    end            
                    
                    % Compute policy function
                    [Vc(iz,ib,ik), pol] = maxmat(V_c_temp);
                    p.ib(iz,ib,ik) = pol(1);
                    p.ik_c(iz,ib,ik) = pol(2);
                    p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                    p.k_c(iz,ib,ik) = grid.k(p.ik_c(iz,ib,ik));
                    p.g_c(iz,ib,ik) = G(pol(1), pol(2));
                                        
                    V1(iz,ib,ik) = max(Vc(iz,ib,ik),Vd(iz,ik));
                    default_states(iz,ib,ik) = 1*(Vd(iz,ik)>Vc(iz,ib,ik));
                    
                    if default_states(iz,ib,ik)==1
                        p.ib(iz,ib,ik) = index0;
                        p.b(iz,ib,ik) = grid.b(index0);
                        p.k(iz,ib,ik) = p.k_d(iz,ik);
                        p.ik(iz,ib,ik) = p.ik_d(iz,ik);
                        p.g(iz,ib,ik) = p.g_d(iz,ik);
                    else
                        p.k(iz,ib,ik) = p.k_c(iz,ib,ik);
                        p.ik(iz,ib,ik) = p.ik_c(iz,ib,ik);
                        p.g(iz,ib,ik) = p.g_c(iz,ib,ik);
                    end
                end
            end
        end
        dist = max(max(max(abs(V1-V))));
        V = V1;
    
        % Update q
        q1 = ones(N.z,N.b,N.k);
        for ik = 1:N.k
            repay = 1-default_states(:,:,ik);
            q1(:,:,ik) = P*repay.*(lambda + (1-lambda)*(coupon+q(:,:,ik)));
        end
        q = q1./(1+r);
        
    end
    save('Policies\GR_ro_1.mat')
    
%% Full commitment beta baseline
clear;close all;clc;
    % Calibration 
    r = 0.01; 
    beta = 0.86; 
    
    sigma = 2;  
    psi = 1/2.22; 
    omega = 0.3; 

    alpha = 0.36;  

    rho = 0.95; 
    sig = 0.01; 

    tau = 0.18;
    delta = 0.025; 
    lambda = 0.05; 
    coupon=0.03;
    zeta = 1.22; 

    % State space
    % Grid points
    N.z = 9;
    N.k = 20;
    N.b = 20;

    % Debt Grid
    bmin = 0;
    bmax = 2.2;
    grid.b = linspace(bmin, bmax, N.b); 
    index0 = 1;

    % Capital Grid
    kmin = 0.15;
    kmax = 3; 
    grid.k = linspace(kmin, kmax, N.k);
    % TFP
    w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
    [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
    grid.z = exp(grid.z)';
    z_mean_index = (N.z+1)/2;
    
    % Initial grids
    V = ones(N.z, N.b, N.k);
    V1 = ones(N.z, N.b, N.k);
    q = (lambda+(1-lambda)*coupon)/(r+lambda);

    
    % VFI loop
    dist = 'Start';
    iter = 0;
    tol = 1e-6;
    while dist>tol        
        iter = iter+1;
        display(num2str([iter dist]))
        
        % Value of staying in debt contract 
        for ib = 1:N.b
            for ik = 1:N.k
                for iz = 1:N.z
                    b = grid.b(ib);
                    l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                    C = y_c/(1+tau);
                    nfa = q.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                    i = grid.k - (1-delta)*grid.k(ik);
                    G = tau*C - i + nfa;
                    neg2 = G<=0  ;
                    uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16; 
             
                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,:,:), N.b, N.k);
                        V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                    end            
                    
                    % Compute policy function
                    [V1(iz,ib,ik), pol] = maxmat(V_c_temp);
                    p.ib(iz,ib,ik) = pol(1);
                    p.ik(iz,ib,ik) = pol(2);
                    p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                    p.k(iz,ib,ik) = grid.k(p.ik(iz,ib,ik));
                    p.g(iz,ib,ik) = G(pol(1),pol(2));
                end
            end
        end
        dist = max(max(max(abs(V1-V))));
        V = V1;
    end
    save('Policies\full_commit_beta_86.mat')
    
    
%% Full commitment beta 0.99
clear;close all;clc;
    % Calibration 
    r = 0.01; 
    beta = 0.99; 
    
    sigma = 2;  
    psi = 1/2.22; 
    omega = 0.3; 

    alpha = 0.36;  

    rho = 0.95; 
    sig = 0.01; 

    tau = 0.18;
    delta = 0.025; 
    lambda = 0.05; 
    coupon=0.03;
    zeta = 1.22; 


    % State space
    % Grid points
    N.z = 9;
    N.k = 20;
    N.b = 20;

    % Debt Grid
    bmin = 0;
    bmax = 2.2;
    grid.b = linspace(bmin, bmax, N.b); 
    index0 = 1;

    % Capital Grid
    kmin = 0.15;
    kmax = 3; 
    grid.k = linspace(kmin, kmax, N.k);
    % TFP
    w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
    [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
    grid.z = exp(grid.z)';
    z_mean_index = (N.z+1)/2;
    
    % Initial grids
    V = ones(N.z, N.b, N.k);
    V1 = ones(N.z, N.b, N.k);
    q = (lambda+(1-lambda)*coupon)/(r+lambda);

    
    % VFI loop
    dist = 'Start';
    iter = 0;
    tol = 1e-6;
    while dist>tol        
        iter = iter+1;
        display(num2str([iter dist]))
        
        % Value of staying in debt contract 
        for ib = 1:N.b
            for ik = 1:N.k
                for iz = 1:N.z
                    b = grid.b(ib);

                    l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                    C = y_c/(1+tau);
                    nfa = q.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                    i = grid.k - (1-delta)*grid.k(ik);
                    
                    G = tau*C - i + nfa;
                    neg2 = G<=0  ;

                    uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16; 
             
                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,:,:), N.b, N.k);
                        V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                    end            
                    
                    % Compute policy function
                    [V1(iz,ib,ik), pol] = maxmat(V_c_temp);
                    p.ib(iz,ib,ik) = pol(1);
                    p.ik(iz,ib,ik) = pol(2);
                    p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                    p.k(iz,ib,ik) = grid.k(p.ik(iz,ib,ik));
                    p.g(iz,ib,ik) = G(pol(1),pol(2));
                end
            end
        end
        dist = max(max(max(abs(V1-V))));
        V = V1;
    end
    save('Policies\full_commit_beta_99.mat')
    
%% Debt ceiling
clear;close all;clc;

    ceiling_grid = [0.3 0.4 0.5];
    for loop_iter=1:numel(ceiling_grid)

        clearvars -except ceiling_grid loop_iter  
        dc = ceiling_grid(loop_iter);   
        % Calibration 
        r = 0.01; 
        beta = 0.86; 
        theta = 0.1; 
        phi = 0.86; 

        sigma = 2; 
        psi = 1/2.22; 
        omega = 0.3; 

        alpha = 0.36; 

        rho = 0.95; 
        sig = 0.01; 

        tau = 0.18; 
        delta = 0.025; 
        lambda = 0.05; 
        coupon=0.03; 
        zeta = 1.22; 

        % State space
        % Grid points
        N.z = 9;
        N.k = 20;
        N.b = 20;

        % Debt Grid
        bmin = 0;
        bmax = 1.1;
        grid.b = linspace(bmin, bmax, N.b); 
        index0 = 1;

        % Capital Grid
        kmin = 0.15;
        kmax = 3; 
        grid.k = linspace(kmin, kmax, N.k);

        % TFP
        w = sig/sqrt(1-rho^2)*(0.5-rho/4)+sig*(0.5+rho/4);
        [grid.z,P] = tauchenhussey(N.z,0,rho,sig,w); 
        grid.z = exp(grid.z)';
        z_mean_index = (N.z+1)/2;


        % Initial grids
        Vd = ones(N.z,N.k);
        Vc = ones(N.z, N.b, N.k);
        V = Vc;
        V1 = Vc;
        default_states = zeros(N.z, N.b, N.k); 
        q =  ones(N.z, N.b, N.k)*(lambda+(1-lambda)*coupon)/(r+lambda);


        % VFI loop
        dist = 'Start';
        iter = 0;
        tol = 1e-6;
        maxit = 1500;
        while dist>tol & iter < maxit        
            iter = iter+1;
            display(num2str([iter dist]))

            % Value of default 
            for ik = 1:N.k
                for iz = 1:N.z
                    k_d = grid.k(ik);
                    z_d = min(grid.z(iz), phi*mean(grid.z));
                    l_d = (((1-alpha)*z_d*zeta*k_d^alpha)/((1+tau)))^(1/(psi+alpha));
                    y_d = zeta*z_d*k_d^alpha*l_d^(1-alpha);
                    C_d =  y_d/(1+tau);
                    i_d =  grid.k - (1-delta)*k_d;
                    G_d = tau*C_d - i_d;
                    neg = G_d<=0;
                    u_d = (C_d -l_d^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                    V_d_temp = (1-omega)*u_d +  omega*G_d.^(1-sigma)/(1-sigma) -neg*10^16;

                    % Compute expectations
                    for iznext = 1:N.z
                        aux = reshape(V(iznext,index0,:), 1, N.k);
                        V_d_temp = V_d_temp + beta*P(iz,iznext).*(theta*aux + (1-theta)*Vd(iznext,:)); 
                    end
                    [Vd(iz,ik), ik_d_opt] = max(V_d_temp);
                    p.ik_d(iz,ik) = ik_d_opt;
                    p.k_d(iz,ik) = grid.k(ik_d_opt);
                    p.g_d(iz,ik) = G_d(ik_d_opt);
                end
            end

            % Value of staying in debt contract 
            for ib = 1:N.b
                for ik = 1:N.k
                    for iz = 1:N.z
                        b = grid.b(ib);
                        qq = q(iz,:,:);
                        qq = reshape(qq, [N.b, N.k]);
                        l_c = (((1-alpha)*grid.z(iz)*zeta*grid.k(ik)^alpha)/((1+tau)))^(1/(psi+alpha));
                        y_c = zeta*grid.z(iz)*(grid.k(ik)^alpha)*l_c^(1-alpha);
                        C = y_c/(1+tau);
                        nfa = qq.*(grid.b'-(1-lambda)*b) - (lambda + (1-lambda)*coupon)*b;
                        i = grid.k - (1-delta)*grid.k(ik);
                        G = tau*C - i + nfa;
                        neg2 = G<=0;
                        constraint_c = grid.b'.*ones(N.b, N.k) > dc*y_c;
                        uc = (C -l_c^(1+psi)/(1+psi)).^(1-sigma)/(1-sigma);
                        V_c_temp = (1-omega)*uc + omega*G.^(1-sigma)/(1-sigma) - neg2*10^16 - constraint_c*10^16; 

                        % Compute expectations
                        for iznext = 1:N.z
                            aux = reshape(V(iznext,:,:), N.b, N.k);
                            V_c_temp = V_c_temp + beta*P(iz,iznext)*aux;
                        end            

                        % Compute policy function
                        [Vc(iz,ib,ik), pol] = maxmat(V_c_temp);
                        p.ib(iz,ib,ik) = pol(1);
                        p.ik_c(iz,ib,ik) = pol(2);
                        p.b(iz,ib,ik) = grid.b(p.ib(iz,ib,ik));
                        p.k_c(iz,ib,ik) = grid.k(p.ik_c(iz,ib,ik));
                        p.g_c(iz,ib,ik) = G(pol(1), pol(2));

                        V1(iz,ib,ik) = max(Vc(iz,ib,ik),Vd(iz,ik));
                        default_states(iz,ib,ik) = 1*(Vd(iz,ik)>Vc(iz,ib,ik));

                        if default_states(iz,ib,ik)==1
                            p.ib(iz,ib,ik) = index0;
                            p.b(iz,ib,ik) = grid.b(index0);
                            p.k(iz,ib,ik) = p.k_d(iz,ik);
                            p.ik(iz,ib,ik) = p.ik_d(iz,ik);
                            p.g(iz,ib,ik) = p.g_d(iz,ik);
                        else
                            p.k(iz,ib,ik) = p.k_c(iz,ib,ik);
                            p.ik(iz,ib,ik) = p.ik_c(iz,ib,ik);
                            p.g(iz,ib,ik) = p.g_c(iz,ib,ik);
                        end
                    end
                end
            end
            dist = max(max(max(abs(V1-V))));
            V = V1;

            % Update q
            q1 = ones(N.z,N.b,N.k);
            for ik = 1:N.k
                repay = 1-default_states(:,:,ik);
                q1(:,:,ik) = P*repay.*(lambda + (1-lambda)*(coupon+q(:,:,ik)));
            end
            q = q1./(1+r);

        end

        str = num2str(dc);
        str = strcat('Policies\DC_', str,'.mat');
        save(str)
    end