classdef RelativeDynamics < handle
    % RELATIVEDYNAMICS Spacecraft Relative Motion Dynamics (Stateful)
    %
    %   This class maintains the state of the chaser internally.
    %   Inherits from 'handle' to allow state updates within methods.
    
    properties
        % Physical Parameters
        J_c      % (3x3) Inertia matrix of the chaser [kg*m^2]
        inv_J_c  % (3x3) Inverse inertia matrix
        m_c      % (1x1) Mass of the chaser [kg]
        MU       % (1x1) Gravitational parameter [m^3/s^2]
        
        % System State (Managed internally)
        state    % (12x1) [sigma; omega; rho; vel]
        time     % Current simulation time [s]
        dt       % Time step
    end
    
    properties (Dependent)
        % Helper properties to access state components easily
        sigma    % (3x1) MRP Attitude
        omega    % (3x1) Angular Velocity
        rho      % (3x1) Relative Position
        vel      % (3x1) Relative Velocity
    end
    
    methods
        % Constructor
        function obj = RelativeDynamics(J_c, m_c, initial_state, dt)
            % J_c: Inertia matrix (3x3)
            % m_c: Mass (scalar)
            % initial_state: (12x1) Initial state vector
            % dt: Dynamics time step
            
            obj.J_c = J_c;
            obj.inv_J_c = inv(J_c);
            obj.m_c = m_c;
            obj.MU = 398600.4418 * 1e9; % Default Earth (m^3/s^2) if not provided
            obj.dt = dt;
            
            % Initialize State
            if nargin < 3 || isempty(initial_state)
                obj.state = zeros(12, 1);
            else
                obj.state = initial_state(:); % Ensure column vector
            end
            obj.time = 0.0;
        end
        
        % Time Update (Main Method)
        function step(obj, u_ctrl, u_dist, target_state)
            % STEP Performs RK4 integration and updates internal state
            %
            % Usage: obj.step(dt, u_ctrl, u_dist, target_state)
            
            t_curr = obj.time;
            x_curr = obj.state;
            
            % RK4 Integration
            k1 = obj.dynamics(           t_curr,                 x_curr,    u_ctrl, u_dist, target_state);
            k2 = obj.dynamics(t_curr + obj.dt/2,   x_curr + obj.dt/2*k1,    u_ctrl, u_dist, target_state);
            k3 = obj.dynamics(t_curr + obj.dt/2,   x_curr + obj.dt/2*k2,    u_ctrl, u_dist, target_state);
            k4 = obj.dynamics(  t_curr + obj.dt,     x_curr + obj.dt*k3,    u_ctrl, u_dist, target_state);
            
            % Update State
            obj.state = x_curr + (obj.dt/6) * (k1 + 2*k2 + 2*k3 + k4);
            
            % Update Time
            obj.time = t_curr + obj.dt;
        end

        % Dependent Property Getters
        function val = get.sigma(obj)
            val = obj.state(1:3);
        end
        function val = get.omega(obj)
            val = obj.state(4:6);
        end
        function val = get.rho(obj)
            val = obj.state(7:9);
        end
        function val = get.vel(obj)
            val = obj.state(10:12);
        end
        
        % Math Helpers
        function S = skew(~, x)
            S = [0,    -x(3),  x(2);
                 x(3),  0,    -x(1);
                -x(2),  x(1),  0];
        end
        
        function G = get_G_matrix(obj, s)
            s_sq = s' * s;
            G = 0.25 * ((1 - s_sq)*eye(3) + 2*obj.skew(s) + 2*(s*s'));
        end
        
        function R = get_Rt_c(obj, s)
            s_sq = s' * s;
            Omega_s = obj.skew(s);
            denom = (1 + s_sq)^2;
            R = eye(3) - (4*(1-s_sq)/denom)*Omega_s + (8/denom)*(Omega_s*Omega_s);
        end
        
        % Dynamics (Internal Calculation)
        function dstate = dynamics(obj, t, x, u_ctrl, u_dist, target_state)
            % Calculates dx/dt given a specific state x.
            % Note: We pass 'x' explicitly to support RK4 intermediate steps.
            
            % Unpacking
            s_curr = x(1:3);
            w_curr = x(4:6);
            r_curr = x(7:9);
            v_curr = x(10:12);
            
            tau = u_ctrl.tau;
            tau_d = u_dist.tau_d;
            f_ctrl = u_ctrl.f;
            f_d = u_dist.f_d;
            
            w_t = target_state.w_t;
            dw_t = target_state.dw_t;
            r_t = target_state.r_t;
            dv_t = target_state.dv_t;
            
            % Kinematics and Dynamics Setup
            R_tc = obj.get_Rt_c(s_curr);
            Rw_t = R_tc * w_t;
            w_c = w_curr + Rw_t; % Chaser angular velocity
            
            Omega_wc = obj.skew(w_c);
            Omega_Rwt = obj.skew(Rw_t);
            
            % dsigma
            dsigma = obj.get_G_matrix(s_curr) * w_curr;
            
            % domega
            C1 = -obj.J_c*Omega_Rwt - Omega_Rwt*obj.J_c + obj.skew(obj.J_c*w_c);
            D1 = -Omega_Rwt*(obj.J_c*Rw_t) - obj.J_c*(R_tc*dw_t);
            
            domega = obj.inv_J_c * ((C1 * w_curr) + D1 + tau + tau_d);
            
            % drho
            drho = v_curr - (Omega_wc * r_curr);
            
            % dv
            r_c = r_curr + (R_tc * r_t); % Chaser position vector
            r_c_norm = norm(r_c);
            
            C2 = -Omega_wc;
            D2 = -(obj.MU / r_c_norm^3) * r_c - (R_tc * dv_t);
            
            dv = ( (obj.m_c * (C2 * v_curr)) + (obj.m_c * D2) + f_ctrl + f_d ) / obj.m_c;
            
            dstate = [dsigma; domega; drho; dv];
        end
        
    end
end