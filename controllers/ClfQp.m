classdef ClfQp < handle
    % Backstepping-based Control Lyapunov Function Quadratic Program (CLF-QP) controller
    %
    % This class implements a CLF-QP controller using a backstepping design for
    % a relative dynamics model (provided via the RelativeDynamics object).
    % The controller computes a control command by first evaluating Lyapunov
    % functions for both position and attitude subsystems (V1 and V2) and then
    % solving a quadratic program to enforce CLF constraints.
    properties
        ref_vel
        prev_ref_vel
        ref_omg
        prev_ref_omg

        rhoV
        velV
        sigV
        omgV

        V1
        V2

        gamma_rho
        gamma_vel
        gamma_sig
        gamma_omg

        torque_lb
        torque_ub
        force_lb
        force_ub

        qp_option
        
        p_weight_f
        p_weight_m

        a_h
        delta_h

        RD
    end

    methods
        function obj = ClfQp(Cfg, relativeDynamics)
            obj.ref_vel = 0;
            obj.prev_ref_vel = NaN;
            obj.ref_omg = 0;
            obj.prev_ref_omg = NaN;

            obj.rhoV = NaN;   % 0.5*(rho_n'*rho_n)
            obj.velV = NaN;   % 0.5*(err_vel'*err_vel) where err_vel = vel - ref_vel
            obj.sigV = NaN;   % 0.5*(sig'*sig)
            obj.omgV = NaN;   % 0.5*(err_omg'*err_omg) where err_omg = omg - ref_omg

            obj.V1 = struct('rel_pos', 0,...
                            'rel_att', 0);          % Backstepping x1
            obj.V2 = struct('rel_pos_vel', 0,...
                            'rel_att_omg', 0);      % Backstepping x2

            obj.gamma_rho  = Cfg.gamma_rho;
            obj.gamma_vel  = Cfg.gamma_vel;
            obj.gamma_sig  = Cfg.gamma_sig;
            obj.gamma_omg  = Cfg.gamma_omg;

            obj.torque_lb = [Cfg.torque_lb(:); 0];
            obj.torque_ub = [Cfg.torque_ub(:); Inf];
            obj.force_lb = [Cfg.force_lb(:); 0];
            obj.force_ub = [Cfg.force_ub(:); Inf];

            obj.qp_option = Cfg.qp_option;

            obj.a_h = Cfg.a_h;
            obj.delta_h = Cfg.delta_h;

            obj.p_weight_f = 1e3;
            obj.p_weight_m = 1e4;

            obj.RD = relativeDynamics;
        end

        function u_ctrl = command(obj)
            obj.ref_vel_cal();
            obj.ref_omg_cal();

            u_ctrl = struct('f', obj.command_force(),...
                            'tau', obj.command_torque());
        end

        function lyapunov_cal(obj)
            sigma = obj.RD.state(1:3);
            omega = obj.RD.state(4:6);
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);

            r_d_t = [obj.delta_h; 0; 0];
            rho_n = rho - obj.RD.R_tc * r_d_t;

            err_vel = vel - obj.ref_vel;
            err_omg = omega - obj.ref_omg;
            obj.rhoV = 0.5 * (rho_n' * rho_n);
            obj.velV = 0.5 * (err_vel' * err_vel);
            obj.sigV = 0.5 * (sigma' * sigma);
            obj.omgV = 0.5 * (err_omg' * err_omg);
            obj.V1.rel_pos = obj.rhoV;
            obj.V1.rel_att = obj.sigV;
            obj.V2.rel_pos_vel = obj.rhoV + obj.velV;
            obj.V2.rel_att_omg = obj.sigV + obj.omgV;
        end

        function ref_vel_cal(obj)
            rho = obj.RD.state(7:9);
            
            r_d_t = [obj.delta_h; 0; 0];
            rho_d_c = obj.RD.R_tc * r_d_t;
            rho_n = rho - rho_d_c;
            
            w_t = obj.RD.Target.stateECI(10:12);
            Rw_t = obj.RD.R_tc * w_t;
            
            v_target_motion = cross(Rw_t, rho_d_c); 
            
            % Lie Derivatives for V_rho = 0.5 * rho_n' * rho_n
            LfV = -rho_n' * v_target_motion; 
            LgV = rho_n'; % (1x3)
            obj.rhoV = 0.5 * (rho_n' * rho_n);

            H = eye(3);
            f = zeros(3,1);

            A = LgV;
            b = -obj.gamma_rho * obj.rhoV - LfV;

            % Solve
            [x_sol, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);

            if exitflag > 0
                obj.ref_vel = x_sol;
            else
                % Fallback
                error('QP solver failed to find a solution for ref_vel');
            end
        end

        function ref_omg_cal(obj)
            sigma = obj.RD.state(1:3);
            G = obj.RD.get_G_matrix(sigma);

            % Lie Derivatives for V_sig = 0.5 * sigma' * sigma
            LfV = 0;
            LgV = sigma' * G;
            obj.sigV = 0.5 * (sigma' * sigma);

            H = eye(3);
            f = zeros(3,1);

            A = LgV;
            b = -obj.gamma_sig * obj.sigV - LfV;

            % Solve
            [x_sol, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);

            if exitflag > 0
                obj.ref_omg = x_sol;
            else
                % Fallback
                error('QP solver failed to find a solution for ref_omg');
            end
        end

        function input_F = command_force(obj)
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            
            r_d_t = [obj.delta_h; 0; 0];
            rho_d_c = obj.RD.R_tc * r_d_t;
            rho_n = rho - rho_d_c;
            
            if isnan(obj.prev_ref_vel(1))
                obj.prev_ref_vel = obj.ref_vel;
            end
            dv_r = (obj.ref_vel - obj.prev_ref_vel) / obj.RD.dt;
            obj.prev_ref_vel = obj.ref_vel;
            err_vel = vel - obj.ref_vel;
            
            w_t = obj.RD.Target.stateECI(10:12);
            Rw_t = obj.RD.R_tc * w_t;
            w_c = obj.RD.state(4:6) + Rw_t; 
            Omega_wc = obj.RD.skew(w_c);
            v_target_motion = cross(Rw_t, rho_d_c);
            
            % Lie Derivative for V_vel = V_rho + 0.5 * err_vel' * err_vel
            LfV_rho = rho_n' * (vel - v_target_motion);
            LfV_vel = err_vel' * (-Omega_wc*vel + obj.RD.gravitational_force() - obj.RD.R_tc * obj.RD.Target.gravitational_force() - dv_r);
            LfV = LfV_rho + LfV_vel;
            LgV = err_vel' / obj.RD.m_c;
            obj.velV = 0.5 * (err_vel' * err_vel);
            obj.V2.rel_pos_vel = obj.rhoV + obj.velV;
            
            H = blkdiag(eye(3), 1e-4);
            f = zeros(4,1);
            f(4) = obj.p_weight_f;
            
            A = [LgV, -1]; 
            b = -obj.gamma_vel * obj.V2.rel_pos_vel - LfV;
            
            [sol, ~, exitflag] = quadprog(H, f, A, b, [], [], obj.force_lb, obj.force_ub, [], obj.qp_option);
            
            if exitflag == -2
                warning('Exit flag -2 occured!')
            elseif isempty(sol)
                error('QP solver failed to find a solution for command_force even with slack variable');
            end
            
            input_F = sol(1:3);
        end

        function input_M = command_torque(obj)
            sigma = obj.RD.state(1:3);
            omega = obj.RD.state(4:6);
            
            if isnan(obj.prev_ref_omg(1)) 
                obj.prev_ref_omg = obj.ref_omg;
            end
            domg_r = (obj.ref_omg - obj.prev_ref_omg) / obj.RD.dt;
            obj.prev_ref_omg = obj.ref_omg;
            err_omg = omega - obj.ref_omg;
            G = obj.RD.get_G_matrix(sigma);
            C1 = obj.RD.get_C1();
            D1 = obj.RD.get_D1();
            
            % Lie Derivatives for V_omg = V_sig + 0.5 * err_omg' * err_omg
            LfV = sigma' * G * omega + err_omg' * (obj.RD.J_c\(C1*omega + D1) - domg_r);
            LgV = err_omg' / obj.RD.J_c;
            obj.omgV = 0.5 * (err_omg' * err_omg);
            obj.V2.rel_att_omg = obj.sigV + obj.omgV;
            
            H = blkdiag(eye(3), 1e-4);
            f = zeros(4,1);
            f(4) = obj.p_weight_m;

            A = [LgV, -1];
            b = -obj.gamma_omg * obj.V2.rel_att_omg - LfV;
            
            [sol , ~, exitflag] = quadprog(H, f, A, b, [], [], obj.torque_lb, obj.torque_ub, [], obj.qp_option);
            
            if exitflag == -2
                warning('Exit flag -2 occured!')
            elseif isempty(sol)
                error('QP solver failed to find a solution for command_force even with slack variable');
            end
            
            input_M = sol(1:3);
        end

        function h = barrier_value(obj)
            r_t = (obj.RD.R_tc')*obj.RD.state(7:9);
            h = obj.a_h*(r_t(1) - obj.delta_h)^3 - r_t(2)^2 - r_t(3)^2;
        end
    end
end