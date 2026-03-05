classdef ClfCbfQp < ClfQp
    % Backstepping-based Control Lyapunov Function and Control Barrier Function Quadratic Program (CLF-CBF-QP) controller
    properties

        alpha_rho
        alpha_vel
        alpha_sig
        alpha_omg

        h_rho
        h_vel
        h_sig
        h_omg

        vel_nom
        vel_safe
        f_nom
        f_safe

        partial
    end
    methods
        function obj = ClfCbfQp(ControlCfg, relativeDynamics)
            obj@ClfQp(ControlCfg, relativeDynamics);

            obj.alpha_rho = ControlCfg.alpha_rho;
            obj.alpha_vel = ControlCfg.alpha_vel;
            obj.alpha_sig = ControlCfg.alpha_sig;
            obj.alpha_omg = ControlCfg.alpha_omg;

            obj.h_rho = NaN;
            obj.h_vel = NaN;
            obj.h_sig = NaN;
            obj.h_omg = NaN;

            obj.vel_nom = NaN;
            obj.f_nom = NaN;

            obj.partial = struct('dh_drho', [], 'd2h_drho2', []);
        end


        function u_ctrl = command(obj)
            obj.partial_rho_t();
            obj.ref_vel_cal();
            obj.ref_omg_cal();
            u_ctrl = struct('f', obj.command_force(),...
                            'tau', obj.command_torque());
        end
        function ref_vel_cal(obj)
            ref_vel_cal@ClfQp(obj);
            obj.vel_nom = obj.ref_vel;
            obj.h_rho = obj.barrier_value();

            H = eye(3);
            f = -obj.vel_nom; 

            Lgh = obj.partial.dh_drho * obj.RD.R_tc';
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);
            
            Lfh = -obj.partial.dh_drho * obj.RD.R_tc' * obj.RD.skew(obj.RD.R_tc*w_t) * rho;
            
            A = -Lgh;
            b = Lfh + obj.alpha_rho * obj.h_rho;
            
            % Solve
            [v_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);
            if exitflag ~= 1
                warning('Velocity CBF-QP did not converge, using nominal control');
                v_safe = obj.vel_nom;
            end
            obj.ref_vel = v_safe;
        end
        
        % function ref_omg_cal(obj)
        %     ref_omg_cal@ClfQp(obj);
        % end

        function input_F = command_force(obj)
            obj.f_nom = command_force@ClfQp(obj);
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            
            h = obj.barrier_value();
            dot_h = obj.partial.dh_drho * obj.RD.R_tc' * (vel - obj.RD.skew(obj.RD.R_tc*w_t)*rho);

            rho_t = obj.RD.R_tc' * rho;
            D2 = obj.RD.gravitational_force() - obj.RD.R_tc * obj.RD.Target.gravitational_force();

            dot_rho_t = obj.RD.R_tc' * (vel - obj.RD.skew(obj.RD.R_tc*w_t)*rho);
            
            Lfh = dot_rho_t' * obj.partial.d2h_drho2 * dot_rho_t...
                + obj.partial.dh_drho * (-2 * obj.RD.skew(w_t) * dot_rho_t - obj.RD.skew(w_t)^2 * rho_t + obj.RD.R_tc' * D2)...
                + obj.alpha_rho * dot_h;
            Lgh = obj.partial.dh_drho * obj.RD.R_tc' ./ obj.RD.m_c;
            h_2 = dot_h + obj.alpha_rho * h;
            
            H = eye(3);
            f = -obj.f_nom;

            A = -Lgh;
            b = Lfh + obj.alpha_vel * h_2;
            
            % Solve
            [force_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], obj.force_lb, obj.force_ub, [], obj.qp_option);
            if exitflag ~= 1
                warning('Force CBF-QP did not converge, using nominal control');
                force_safe = obj.f_nom;
            end
            input_F = force_safe;
            obj.f_safe = input_F;
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

            % Lie Derivatives for V_omg = V_sig + 0.5 * err_omg' * err_omg
            LfV = sigma' * G * omega + err_omg' * (obj.RD.J_c\(C1*omega + D1) - domg_r);
            LgV = err_omg' / obj.RD.J_c;
            obj.omgV = 0.5 * (err_omg' * err_omg);
            obj.V2.rel_att_omg = obj.sigV + obj.omgV;

            H = eye(3);
            f = zeros(3,1);

            A = LgV;
            b = -obj.gamma_omg * obj.V2.rel_att_omg - LfV;

            [input_M , ~, exitflag] = quadprog(H, f, A, b, [], [], obj.torque_lb, obj.torque_ub, [], obj.qp_option);
            
            if exitflag == -2
                % Fallback
                error('QP solver failed to find a solution for command_torque');
            end
        end

        function partial_rho_t(obj)
            rho = obj.RD.state(7:9);
            rho_t = obj.RD.R_tc' * rho;
            
            grad_rho_t = [3*obj.a_h*(rho_t(1) - obj.delta_h)^2,...
                       -2*rho_t(2),...
                       -2*rho_t(3)];
            hessian_rho_t = [6*obj.a_h*(rho_t(1) - obj.delta_h), 0, 0;...
                         0, -2, 0;...
                         0, 0, -2];
            obj.partial.dh_drho = grad_rho_t;
            obj.partial.d2h_drho2 = hessian_rho_t;
        end
    end
end