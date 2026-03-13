classdef Backstepping < handle
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

        force_slack
        torque_slack

        qp_option

        a_h
        delta_h

        RD  RelativeDynamics
    end

    methods
        function obj = Backstepping(ControlCfg, relativeDynamics)
            obj.ref_vel = 0;
            obj.prev_ref_vel = NaN;
            obj.ref_omg = 0;
            obj.prev_ref_omg = NaN;

            obj.rhoV = 0;   % 0.5*(rho'*rho)
            obj.velV = 0;   % 0.5*(err_vel'*err_vel) where err_vel = vel - ref_vel
            obj.sigV = 0;   % 0.5*(sig'*sig)
            obj.omgV = 0;   % 0.5*(err_omg'*err_omg) where err_omg = omg - ref_omg

            obj.V1 = struct('rel_pos', 0,...
                'rel_att', 0);          % Backstepping x1
            obj.V2 = struct('rel_pos_vel', 0,...
                'rel_att_omg', 0);      % Backstepping x2

            obj.gamma_rho   = ControlCfg.gamma_rho;
            obj.gamma_vel   = ControlCfg.gamma_vel;
            obj.gamma_sig   = ControlCfg.gamma_sig;
            obj.gamma_omg   = ControlCfg.gamma_omg;

            obj.torque_lb = ControlCfg.torque_lb;
            obj.torque_ub = ControlCfg.torque_ub;
            obj.force_lb = ControlCfg.force_lb;
            obj.force_ub = ControlCfg.force_ub;

            obj.force_slack = ControlCfg.force_slack;
            obj.torque_slack = ControlCfg.torque_slack;

            obj.qp_option = optimoptions('quadprog', 'Display', 'off');

            obj.a_h = ControlCfg.a_h;
            obj.delta_h = ControlCfg.delta_h;

            obj.RD = relativeDynamics;
        end

        function u_ctrl = command(obj)
            obj.lyapunov_cal();
            obj.ref_vel_cal();
            obj.ref_omg_cal();
       
            u_ctrl = struct('f', obj.command_force(),...
                            'tau', obj.command_torque());
            % u_ctrl = obj.saturate(u_ctrl);
        end

        function lyapunov_cal(obj)
            sigma = obj.RD.state(1:3);
            omega = obj.RD.state(4:6);
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            err_vel = vel - obj.ref_vel;
            err_omg = omega - obj.ref_omg;

            obj.rhoV = 0.5 * (rho' * rho);
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

            A = LgV;
            b = -obj.gamma_rho * obj.rhoV - LfV;

            obj.ref_vel = pinv(A)*b;
        end

        function ref_omg_cal(obj)
            sigma = obj.RD.state(1:3);
            G = obj.RD.get_G_matrix(sigma);

            LfV = 0;
            LgV = sigma' * G;

            % Constraint: LgV * ref_omg <= -gamma1 * V - LfV
            % No slack varialbe for now
            A = LgV;
            b = -obj.gamma_sig * obj.sigV - LfV;

            obj.ref_omg = pinv(A)*b;
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

            A = LgV;
            b = -obj.gamma_vel * obj.V2.rel_pos_vel - LfV;
            
            input_F = pinv(A)*b;
        end

        function input_M = command_torque(obj)
            sigma = obj.RD.state(1:3);
            omega = obj.RD.state(4:6);
            if isnan(obj.prev_ref_omg)
                obj.prev_ref_omg = obj.ref_omg;
            end
            domg_r = (obj.ref_omg - obj.prev_ref_omg) / obj.RD.dt;
            obj.prev_ref_omg = obj.ref_omg;
            err_omg = omega - obj.ref_omg;

            G = obj.RD.get_G_matrix(sigma);
            C1 = obj.RD.get_C1();
            D1 = obj.RD.get_D1();
            LfV = sigma' * G * omega + err_omg' * (obj.RD.J_c\(C1*omega + D1) - domg_r);
            % LfV = sigma' * G * omega + err_omg' * (obj.RD.J_c\(C1*omega + D1));
            LgV = err_omg' / obj.RD.J_c;

            A = LgV;
            b = -obj.gamma_omg * obj.V2.rel_att_omg - LfV;
            input_M = pinv(A)*b;
        end

        function FM = saturate(obj, FM)
            for i = 1:3
                if FM.f(i) > obj.force_ub(i)
                    FM.f(i) = obj.force_ub(i);
                elseif FM.f(i) < obj.force_lb(i)
                    FM.f(i) = obj.force_lb(i);
                end

                if FM.tau(i) > obj.torque_ub(i)
                    FM.tau(i) = obj.torque_ub(i);
                elseif FM.tau(i) < obj.torque_lb(i)
                    FM.tau(i) = obj.torque_lb(i);
                end
            end
        end

        function h = barrier_value(obj)
            r_t = (obj.RD.R_tc')*obj.RD.state(7:9);
            h = obj.a_h*(r_t(1) - obj.delta_h)^3 - r_t(2)^2 - r_t(3)^2;
        end
    end
end