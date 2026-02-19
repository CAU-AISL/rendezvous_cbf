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
        function obj = ClfQp(ControlCfg, relativeDynamics)
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

            % Solving QP
            % Variables: x = ref_vel (3x1)
            H = eye(3);
            f = zeros(3,1);

            % Constraint: LgV * ref_vel <= -gamma1 * V - LfV
            % No slack varialbe for now
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

            LfV = 0;
            LgV = sigma' * G;

            % Solving QP
            % Variables: x = ref_omg (3x1)
            H = eye(3);
            f = zeros(3,1);

            % Constraint: LgV * ref_omg <= -gamma1 * V - LfV
            % No slack varialbe for now
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

            % H = blkdiag(eye(3), obj.force_slack);
            H = eye(3);
            % f = zeros(4,1);
            f = zeros(3,1);

            % A = [LgV, 1];
            A = LgV;
            b = -obj.gamma_rho_vel * obj.V2.rel_pos_vel - LfV;
            
            [x_sol, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);
            command_force = zeros(3,1);
            if exitflag > 0
                command_force = x_sol(1:3);
            elseif exitflag == -2
                % Fallback
                error('QP solver failed to find a solution for command_force');
            end
            % for i = 1:3
            %     if A(i)*x_sol(i) > b
            %         command_force(i) = obj.force_lb(i);
            %     elseif A(i)*x_sol(i) < b
            %         command_force(i) = obj.force_ub(i);
            %     else
            %         command_force(i) = x_sol(i);
            %     end
            % end
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

            % H = blkdiag(eye(3), obj.torque_slack);
            H = eye(3);
            % f = zeros(4,1);
            f = zeros(3,1);

            % A = [LgV, 1];
            A = LgV;
            b = -obj.gamma_sig_omg * obj.V2.rel_att_omg - LfV;
            
            % [x_sol, ~, exitflag] = quadprog(H, f, A, b, [], [], obj.torque_lb, obj.torque_ub, [], obj.qp_option);
            [x_sol, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);
            command_torque = zeros(3,1);
            if exitflag > 0
                command_torque = x_sol(1:3);
            elseif exitflag == -2
                % Fallback
                error('QP solver failed to find a solution for command_torque');
            end
            % for i = 1:3
            %     if A(i)*x_sol(i) > b
            %         command_torque(i) = obj.torque_lb(i);
            %     elseif A(i)*x_sol(i) < b
            %         command_torque(i) = obj.torque_ub(i);
            %     else
            %         command_torque(i) = x_sol(i);
            %     end
            % end
        end
    end
end