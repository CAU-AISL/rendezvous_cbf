classdef StochasticRelative < RelativeDynamics
    properties
        sig_rho
        sig_vel
        sig_mrp
        sig_omg
        mu_rho
        mu_vel
        mu_mrp
        mu_omg
        
        % LPF 객체 및 백색 노이즈 스펙
        lpf_rho
        lpf_vel
        lpf_mrp
        lpf_omg
        w_sig_rho
        w_sig_vel
        w_sig_mrp
        w_sig_omg
        
    end
    methods
        function obj = StochasticRelative(Cfg, targetSatellite)
            obj@RelativeDynamics(Cfg, targetSatellite);
            
            obj.sig_rho = Cfg.sigma.rho;
            obj.sig_vel = Cfg.sigma.vel;
            obj.sig_mrp = Cfg.sigma.mrp;
            obj.sig_omg = Cfg.sigma.omg;
            obj.mu_rho = Cfg.mu.rho;
            obj.mu_vel = Cfg.mu.vel;
            obj.mu_mrp = Cfg.mu.mrp;
            obj.mu_omg = Cfg.mu.omg;
            
            % 1st order LPF Initialization
            omg_c = 10;
            obj.lpf_rho = LPF(omg_c, Cfg);
            obj.lpf_vel = LPF(omg_c, Cfg);
            obj.lpf_mrp = LPF(omg_c, Cfg);
            obj.lpf_omg = LPF(omg_c, Cfg);
            
            % Distributed amplification of input white noise by 
            % back-calculating the distributed attenuation of the LPF
            alpha = obj.lpf_rho.alpha; 
            scale_factor = sqrt((1 + alpha) / (1 - alpha));
            obj.w_sig_rho = obj.sig_rho * scale_factor;
            obj.w_sig_vel = obj.sig_vel * scale_factor;
            obj.w_sig_mrp = obj.sig_mrp * scale_factor;
            obj.w_sig_omg = obj.sig_omg * scale_factor;
        end
        
        function d_state = dynamics(obj, state, u_ctrl, u_dist)
            d_state = dynamics@RelativeDynamics(obj, state, u_ctrl, u_dist);
            % Gaussian noise
            w_rho = normrnd(obj.mu_rho, obj.w_sig_rho, [3, 1]);
            w_vel = normrnd(obj.mu_vel, obj.w_sig_vel, [3, 1]);
            w_mrp = normrnd(obj.mu_mrp, obj.w_sig_mrp, [3, 1]);
            w_omg = normrnd(obj.mu_omg, obj.w_sig_omg, [3, 1]);

            
            % Colored Gaussian Noise after 1st order LPF
            c_rho = obj.lpf_rho.forward(w_rho);
            c_vel = obj.lpf_vel.forward(w_vel);
            c_mrp = obj.lpf_mrp.forward(w_mrp);
            c_omg = obj.lpf_omg.forward(w_omg);
            
            d_state = d_state + [c_mrp; c_omg; c_rho; c_vel];
        end
    end
end