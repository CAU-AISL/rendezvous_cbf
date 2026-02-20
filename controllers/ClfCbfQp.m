classdef ClfCbfQp < ClfQp
    % Backstepping-based Control Lyapunov Function and Control Barrier Function Quadratic Program (CLF-CBF-QP) controller
    properties
        torque_lb
        torque_ub
        force_lb
        force_ub

        alpha_rho
        alpha_vel
        alpha_sig
        alpha_omg

        h_rho
        h_vel
        h_sig
        h_omg
    end
    methods
        function obj = ClfCbfQp(ControlCfg, relativeDynamics)
            obj@ClfQp(ControlCfg, relativeDynamics);

            obj.torque_lb = ControlCfg.torque_lb;
            obj.torque_ub = ControlCfg.torque_ub;
            obj.force_lb = ControlCfg.force_lb;
            obj.force_ub = ControlCfg.force_ub;

            obj.alpha_rho = ControlCfg.alpha_rho;
            obj.alpha_vel = ControlCfg.alpha_vel;
            obj.alpha_sig = ControlCfg.alpha_sig;
            obj.alpha_omg = ControlCfg.alpha_omg;

            obj.h_rho = NaN;
            obj.h_vel = NaN;
            obj.h_sig = NaN;
            obj.h_omg = NaN;
        end

        function u_ctrl = command(obj)
            % Cbf constraint cascade or HOCBF
            % Consider input constraints
        end

        function cbf_cal(obj)
            % obj.h_rho
            % obj.h_rho
            % obj.h_rho
        end

        function ref_vel_cal(obj)
            ref_vel_cal@ClfQp(obj);
            % obj.h_rho = ...
        end
        
        function ref_omg_cal(obj)
            ref_omg_cal@ClfQp(obj);
            % obj.h_sig = ...
        end

        function input_F = command_force(obj)
            input_F = command_force@ClfQp(obj);
            % obj.h_vel = ...
        end

        function input_M = command_torque(obj)
            input_M = command_torque@ClfQp(obj);
            % obj.h_omg = ...
        end
    end
end