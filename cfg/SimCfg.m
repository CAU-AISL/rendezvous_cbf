classdef SimCfg < handle
    properties
        % --- Simulation Configuration ---
        dt
        T
        sim_time
        sim_len
        target_init_state
        chaser_init_state

        % --- Control Configuration ---
        gamma_rho
        gamma_vel
        gamma_sig
        gamma_omg
        
        torque_lb
        torque_ub
        force_lb
        force_ub

        force_slack
        torque_slack

        alpha_rho
        alpha_vel
        alpha_sig
        alpha_omg

        eta_rho
        eta_vel

        a_h
        delta_h

        sigma
        mu
        
        qp_option
    end
    methods
        function obj = SimCfg()
            % ==========================================
            % Simulation Configuration
            % ==========================================
            obj.dt  = 0.01;
            obj.T   = 200;
            obj.sim_time = 0 : obj.dt : obj.T;
            obj.sim_len = length(obj.sim_time);
            obj.target_init_state = struct('a', 7702455,...
                                           'e', 0.12,...
                                           'i', deg2rad(30),...
                                           'Omega', deg2rad(0),...
                                           'omega', deg2rad(0),...
                                           'f0', deg2rad(0));
            obj.chaser_init_state =  struct('sigma', [-0.1; 0.12; 0.1],...
                                            'omega', [0.05; -0.03; 0.07],...
                                            'rho', [47.2; -16.6; 38.4],...
                                            'vel', [-0.2; -0.3; -0.1]);

            % ==========================================
            % Control Configuration
            % ==========================================
            obj.gamma_rho = 0.8;
            obj.gamma_vel = 0.08;
            obj.gamma_sig = 3;
            obj.gamma_omg = 0.1;

            obj.torque_lb = [-5; -5; -5];
            obj.torque_ub = [5; 5; 5];
            obj.force_lb = [-20; -20; -20];
            obj.force_ub = [20; 20; 20];

            obj.alpha_rho = 0.8;
            obj.alpha_vel = 0.1;
            obj.alpha_sig = 1;
            obj.alpha_omg = 1;

            obj.eta_rho = 0.8;
            obj.eta_vel = 0.9;

            obj.a_h = 0.1;
            obj.delta_h = 1;

            obj.sigma = struct('rho', 0.01*ones([3, 1]),...
                               'vel', 0.005*ones([3, 1]),...
                               'mrp', deg2rad(0.1)*ones([3, 1]),...
                               'omg', deg2rad(0.03)*ones([3, 1]));
            obj.mu = struct('rho', [0.01; -0.02; 0.03],...
                            'vel', [0.02; 0.01; -0.01],...
                            'mrp', 1e-5*[2; -1; 5],...
                            'omg', 1e-4*[-3; 4; -1]);

            obj.qp_option = optimoptions('quadprog', 'Display', 'off', 'ConstraintTolerance', 1e-5);
        end
    end
end