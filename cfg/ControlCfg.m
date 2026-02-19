classdef ControlCfg < handle
    properties
        gamma_rho
        gamma_rho_vel
        gamma_sig
        gamma_sig_omg
        
        torque_lb
        torque_ub
        force_lb
        force_ub

        force_slack
        torque_slack
    end
    methods
        function obj = ControlCfg()
            obj.gamma_rho = 5;
            obj.gamma_rho_vel = 0.5;
            obj.gamma_sig = 5;
            obj.gamma_sig_omg = 0.2;

            obj.torque_lb = [-5; -5; -5];
            obj.torque_ub = [5; 5; 5];
            obj.force_lb = [-20; -20; -20];
            obj.force_ub = [20; 20; 20];

            obj.force_slack = 1000;
            obj.torque_slack = 1000;
        end
    end
end