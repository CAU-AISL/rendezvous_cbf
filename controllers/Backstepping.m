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
        gamma_rho_vel
        gamma_sig
        gamma_sig_omg
        
        torque_lb
        torque_ub
        force_lb
        force_ub

        force_slack
        torque_slack

        qp_option

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

            obj.gamma_rho      = ControlCfg.gamma_rho;
            obj.gamma_rho_vel  = ControlCfg.gamma_rho_vel;
            obj.gamma_sig      = ControlCfg.gamma_sig;
            obj.gamma_sig_omg  = ControlCfg.gamma_sig_omg;

            obj.torque_lb = ControlCfg.torque_lb;
            obj.torque_ub = ControlCfg.torque_ub;
            obj.force_lb = ControlCfg.force_lb;
            obj.force_ub = ControlCfg.force_ub;

            obj.force_slack = ControlCfg.force_slack;
            obj.torque_slack = ControlCfg.torque_slack;

            obj.qp_option = optimoptions('quadprog', 'Display', 'off');

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

            % Lie Derivatives for V_rho = 0.5 * rho' * rho
            % dot(V) = rho' * (v - S(w)rho) = rho' * v  (since rho'*S*rho = 0)
            LfV = 0;
            LgV = rho'; % (1x3)

            % Constraint: LgV * ref_vel <= -gamma1 * V - LfV
            % No slack varialbe for now
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

        function command_force = command_force(obj)
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            if isnan(obj.prev_ref_vel)
                obj.prev_ref_vel = obj.ref_vel;
            end
            dv_r = (obj.ref_vel - obj.prev_ref_vel) / obj.RD.dt;
            obj.prev_ref_vel = obj.ref_vel;
            err_vel = vel - obj.ref_vel;

            R_tc = obj.RD.get_Rt_c(obj.RD.state(1:3));
            w_t = obj.RD.Target.stateECI(10:12);
            Rw_t = R_tc * w_t;
            w_c = obj.RD.state(4:6) + Rw_t; % Chaser relative angular velocity
            Omega_wc = obj.RD.skew(w_c);
            obj.RD.skew(obj.RD.state(4:6));
            LfV = rho' * vel + err_vel'*(-Omega_wc*vel + obj.RD.gravitational_force() - R_tc * obj.RD.Target.gravitational_force() - dv_r);
            % LfV = rho' * vel + err_vel'*(-Omega_wc*vel + obj.RD.gravitational_force() - R_tc * obj.RD.Target.gravitational_force());
            LgV = err_vel'/obj.RD.m_c;

            % A = [LgV, 1];
            A = LgV;
            b = -obj.gamma_rho_vel * obj.V2.rel_pos_vel - LfV;
            command_force = pinv(A)*b;
        end

        function command_torque = command_torque(obj)
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
            b = -obj.gamma_sig_omg * obj.V2.rel_att_omg - LfV;
            command_torque = pinv(A)*b;
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
    end
end