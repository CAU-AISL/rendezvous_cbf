classdef NominalVel < handle
    % NominalVel Implements CLF-CBF-QP based nominal controller
    % Incorporates coordinate transformation for CBF defined in Target Frame
    
    properties
        gamma   % CLF decay rate
        alpha   % CBF decay rate (Class K function)
        p_relax % High cost for slack variable
        
        qp_options % Quadprog options
        
        Vd_prev % Previous desired velocity command
        dt      % Time step
        V       % Value of CLF
        LfV     % Lie derivative of V along f
        LgV     % Lie derivative of V along g
        B       % Value of CBF
        LfB     % Lie derivative of h along f
        LgB     % Lie derivative of h along g
        
        RelativeChaser  RelativeDynamics
    end
    
    methods
        function obj = NominalVel(control_setting, relativeDynamics)
            obj.gamma = control_setting.gamma;
            obj.alpha = control_setting.alpha;
            obj.p_relax = control_setting.p_relax;
            
            obj.qp_options = control_setting.qp_options;
            
            obj.Vd_prev = NaN;
            obj.dt = relativeDynamics.dt;
            % Initialize properties
            obj.V   = 0;
            obj.LfV = 0;
            obj.LgV = zeros([1, 3]);
            obj.B   = 0;
            obj.LfB = 0;
            obj.LgB = zeros([1, 3]);
            

            obj.RelativeChaser = relativeDynamics;
        end

        function [Force, slack, feas] = command(obj, V_d)
            state = obj.RelativeChaser.state;
            obj.update_clf(state, V_d);
            % obj.update_cbf(state);
            
            quadprog_H = eye(3);
            quadprog_f = zeros([3, 1]);

            A = obj.LgV;
            b = -obj.gamma * obj.V - obj.LfV;
            [output, ~, exitflag, ~] = quadprog(quadprog_H, quadprog_f, A, b, [], [], [], [], [], obj.qp_options);
            
            if exitflag == -2
                Force = zeros(3,1);
                feas = 0;
                disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
            else
                Force = output;
                feas = 1;
            end
            slack = 0;
        end
        
        % ========================================================
        %%%%%%%%%% CLF Functions (V = 0.5 * eV' * eV) %%%%%%%%%%
        % ========================================================
        function update_clf(obj, state, V_d)
            % Compute value of control lyapunov function V
            vel = state(10:12);
            eV = V_d - vel;
            obj.V = 0.5 * (eV' * eV);
            
            % Compute LfV and LgV
            if isnan(obj.Vd_prev)
                obj.Vd_prev = V_d;
            end
            dotVd = (vel - obj.Vd_prev)./obj.dt;
            obj.Vd_prev = V_d;

            S_omega_c = obj.RelativeChaser.skew(obj.RelativeChaser.omg_c);
            obj.LfV = eV' * (S_omega_c * vel + obj.RelativeChaser.Target.gravitational_force() + obj.RelativeChaser.gravitational_force());
            % obj.LfV = eV' * (dotVd + S_omega_c * vel + obj.RelativeChaser.Target.gravitational_force() + obj.RelativeChaser.gravitational_force());
            obj.LgV = -eV'./obj.RelativeChaser.m_c;
        end

        
        % =========================================================
        % CBF Functions (h in Target Frame)
        % =========================================================

        % function update_cbf(obj, state)
        %     % Compute cbf h
        %     R_tc = obj.RelativeChaser.get_Rt_c(state(1:3));
        %     r_t = R_tc' * state(7:9); 
            
        %     obj.B = obj.Alpha * (r_t(1) + obj.Delta)^3 - r_t(2)^2 - r_t(3)^2;
            
        %     % Compute Lfh and Lgh
        %     dh_dxt = 3 * obj.Alpha * (r_t(1) + obj.Delta)^2;
        %     dh_dyt = -2 * r_t(2);
        %     dh_dzt = -2 * r_t(3);
        %     grad_h_rt = [dh_dxt, dh_dyt, dh_dzt];
            
        %     S_omega = obj.RelativeChaser.skew(state(4:6));
        %     S_omega_c = obj.RelativeChaser.skew(obj.RelativeChaser.omg_c);
        %     rho = state(7:9);

        %     obj.LfB = grad_h_rt * (-S_omega * R_tc' * rho - R_tc' * S_omega_c * rho);
        %     obj.LgB = grad_h_rt * R_tc';
        % end
    end
end