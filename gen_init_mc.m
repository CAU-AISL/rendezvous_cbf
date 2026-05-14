clear all; close all; clc;
N = 100;

chaser_init_state = repmat(struct('sigma', zeros(3,1), ...
                                  'omega', zeros(3,1), ...
                                  'rho', zeros(3,1), ...
                                  'vel', zeros(3,1)), N, 1);

count = 1;
while count <= N
    sigma = deg2rad(-5) + rand(3, 1) * deg2rad(10);
    omega = deg2rad(-3) + rand(3, 1) * deg2rad(6);
    
    S_sigma = [0, -sigma(3), sigma(2);
               sigma(3), 0, -sigma(1);
               -sigma(2), sigma(1), 0];
    
    sig_norm2 = sigma' * sigma;
    denom = (1 + sig_norm2)^2;
    
    R_t_c = eye(3) - (4 * (1 - sig_norm2) / denom) * S_sigma + (8 / denom) * (S_sigma * S_sigma);
    
    valid_rho = false;
    while ~valid_rho
    
        u = randn(3, 1);
        u = u / norm(u);
    
        rho_t = 50 * u;
        
        xt = rho_t(1);
        yt = rho_t(2);
        zt = rho_t(3);
    
        if 0.1 * (xt - 1)^3 - yt^2 - zt^2 >= 0
            valid_rho = true;
        end
    end
    
    rho = R_t_c * rho_t;
    vel = -0.1 + 0.2*rand(3, 1); 
    
    chaser_init_state(count).sigma = sigma;
    chaser_init_state(count).omega = omega;
    chaser_init_state(count).rho = rho;
    chaser_init_state(count).vel = vel;
    
    count = count + 1;
end
save(fullfile(pwd, 'MonteCarloInitCondition'), 'chaser_init_state');
