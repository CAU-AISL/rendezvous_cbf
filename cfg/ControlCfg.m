classdef ControlCfg < handle
    properties
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

        a_h
        delta_h

        qp_option
    end

    methods
        function obj = ControlCfg()
            obj.gamma_rho = 0.8;
            obj.gamma_vel = 0.08;
            obj.gamma_sig = 3;
            obj.gamma_omg = 0.1;

            obj.torque_lb = [-5; -5; -5];
            obj.torque_ub = [5; 5; 5];
            obj.force_lb = [-20; -20; -20];
            obj.force_ub = [20; 20; 20];

            % obj.alpha_rho = 0.2;
            % obj.alpha_vel = 0.2;
            obj.alpha_rho = 0.8;
            obj.alpha_vel = 0.1;
            obj.alpha_sig = 1;
            obj.alpha_omg = 1;

            obj.a_h = 0.1;
            obj.delta_h = 1;

            obj.qp_option = optimoptions('quadprog', 'Display', 'off', 'ConstraintTolerance', 1e-5);
        end
    end
end