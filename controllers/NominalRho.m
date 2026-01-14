classdef NominalRho < handle
    % NominalRho Implements CLF-CBF-QP based nominal controller
    % This controller uses backstepping based approach and makes desired relative velocity

    properties
        gamma   % CLF decay rate
        alpha   % CBF decay rate
        p_relax % High cost for slack variable
        u_lim  % Control input limits
        
        qp_options % Quadprog options
        
        V       % Value of CLF
        LfV     % Lie derivative of V along f
        LgV     % Lie derivative of V along g
        b       % Value of CBF
        Lfb     % Lie derivative of b along f
        Lgb     % Lie derivative of h along g
    end

    methods
        function obj = NominalRho(control_setting)
            obj.gamma = control_setting.gamma;
            obj.alpha = control_setting.alpha;
            obj.p_relax = control_setting.p_relax;
            obj.u_lim = control_setting.u_lim;
            
            obj.qp_options = control_setting.qp_options;

            obj.V   = 0;
            obj.LfV = 0;
            obj.LgV = 0;
            obj.b   = 0;
            obj.Lfb = 0;
            obj.Lgb = 0;
        end

        function compute_clf(obj, state)
            rho = state(7:9);
            obj.V = 0.5 * (rho' * rho);
        end

        function compute_LfV(obj, state)
            rho = state(7:9);
            vel = state(10:12);

            obj.LfV = rho' * (vel);
        end

        function compute_LgV(obj, state)
            obj.LgV = state(7:9)';
        end

        function [u_ctrl, debug_info] = compute_control(obj, rel_dyn)
            % compute_control Calculates control inputs using CLF-CBF-QP
            % rel_dyn: Instance of RelativeDynamics class

            x = rel_dyn.state;

            % Get Dynamics: dx = f(x) + g(x)u
            [f, g] = obj.get_affine_dynamics(rel_dyn);

            % 2. CLF Formulation: V = x'Px
            % Objective: Stabilize x to 0.
            % Constraint: LfV + LgV*u + delta <= -gamma*V
            % Rearranged: LgV*u - delta <= -gamma*V - LfV

            V = x' * obj.P * x;
            LfV = x' * (obj.P * f + obj.P' * f);
            LgV = 2 * x' * obj.P * g; % (1x6 vector)

            A_clf = [LgV, -1];
            b_clf = -obj.gamma * V - LfV;

            % 3. CBF Formulation: h(x) >= 0 (Max Velocity Constraint)
            % h(x) = v_max^2 - ||v||^2 >= 0
            % dh/dt = -2 v' * dv/dt = -2 v' (f_v + g_v u)
            % Constraint: dh/dt + alpha*h >= 0
            % Rearranged: -Lgh*u <= alpha*h + Lfh

            v_curr = x(10:12);
            h = obj.v_max^2 - (v_curr' * v_curr);

            % Extract velocity components (indices 10:12)
            f_v = f(10:12);
            g_v = g(10:12, :);

            Lfh = -2 * v_curr' * f_v;
            Lgh = -2 * v_curr' * g_v;

            % Note: A*x <= b.
            % Condition: Lfh + Lgh*u >= -alpha*h
            % -Lgh*u <= alpha*h + Lfh

            A_cbf = [-Lgh, 0]; % Slack variable not used in CBF here (hard constraint)
            b_cbf = obj.alpha * h + Lfh;

            % 4. Solve QP
            % Variables: U = [u_tau (3); u_force (3); delta (1)]
            % Objective: min u'u + p_relax * delta^2

            H = diag([ones(6,1); obj.p_relax]);
            F = zeros(7, 1);

            A = [A_clf; A_cbf];
            b = [b_clf; b_cbf];

            % Input Limits
            lb = [-obj.u_max; 0];
            ub = [obj.u_max; Inf];

            try
                [u_aug, ~, exitflag] = quadprog(H, F, A, b, [], [], lb, ub, [], obj.qp_options);
            catch
                exitflag = -2; % Quadprog error
            end

            if exitflag == 1
                u_val = u_aug(1:6);
                delta_val = u_aug(7);
            else
                % Fallback (e.g. if QP fails, do nothing or passive safety)
                % warning('QP Failed with exitflag %d', exitflag);
                u_val = zeros(6,1);
                delta_val = 0;
            end

            % Output
            u_ctrl.tau = u_val(1:3);
            u_ctrl.f = u_val(4:6);

            debug_info.V = V;
            debug_info.h = h;
            debug_info.delta = delta_val;
            debug_info.exitflag = exitflag;
        end

        function [f, g] = get_affine_dynamics(~, rel_dyn)
            % Numerically extract f and g from rel_dyn: dx = f(x) + g(x)u
            % This ensures consistency with the simulation model.

            x = rel_dyn.state;

            % Zero inputs
            u0.tau = zeros(3,1);
            u0.f = zeros(3,1);
            ud0.tau_d = zeros(3,1);
            ud0.f_d = zeros(3,1);

            % Get drift dynamics f(x)
            f = rel_dyn.dynamics(x, u0, ud0);

            % Get input matrix g(x)
            g = zeros(12, 6);

            % Since dynamics are affine control, we can perturb by 1.0 to get columns
            % Column 1-3: Torque inputs
            for i=1:3
                ut = u0; ut.tau(i) = 1.0;
                dx = rel_dyn.dynamics(x, ut, ud0);
                g(:, i) = dx - f;
            end

            % Column 4-6: Force inputs
            for i=1:3
                uf = u0; uf.f(i) = 1.0;
                dx = rel_dyn.dynamics(x, uf, ud0);
                g(:, 3+i) = dx - f;
            end
        end
    end
end
