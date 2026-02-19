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
        state   % (12x1) [sigma; omega; rho; vel]
        dt      % Time step
        rho_c   % Chaser position in ECI frame
        att_c   % Chaser attitude relative to ECI
        omg_c   % Chaser angular velocity relative to ECI

        Target  SatelliteDynamics
    end
    
    methods
        % Constructor
        function obj = RelativeDynamics(initial_state, dt, targetSatellite)
            % J_c: Inertia matrix (3x3)
            % m_c: Mass (scalar)
            % initial_state: (12x1) Initial state vector
            % dt: Dynamics time step
            
            obj.J_c = [124.4, 22.5, -21.5;...
                        22.5, 163.6, -7;...
                       -21.5, -7, 128.3];
            obj.inv_J_c = eye(3)/obj.J_c;
            obj.m_c = 38.2;
            obj.MU = 3.986004e14;
            obj.dt = dt;
            
            obj.state = zeros([12, 1]);
            obj.state(1:3) = initial_state.sigma;
            obj.state(4:6) = initial_state.omega;
            obj.state(7:9) = initial_state.rho;
            obj.state(10:12) = initial_state.vel;

            obj.rho_c = zeros([3, 1]);
            obj.att_c = zeros([3, 1]);
            obj.omg_c = zeros([3, 1]);

            obj.Target = targetSatellite;
        end
        
        % Time Update (Main Method)
        function step(obj, u_ctrl, u_dist)
            % STEP Performs RK4 integration and updates internal state
            %
            % Usage: obj.step(dt, u_ctrl, u_dist, target_state)
            
            x_curr = obj.state;
            
            % RK4 Integration
            k1 = obj.dynamics(              x_curr, u_ctrl, u_dist);
            k2 = obj.dynamics(x_curr + obj.dt/2*k1, u_ctrl, u_dist);
            k3 = obj.dynamics(x_curr + obj.dt/2*k2, u_ctrl, u_dist);
            k4 = obj.dynamics(  x_curr + obj.dt*k3, u_ctrl, u_dist);
            
            % Update State
            obj.state = x_curr + (obj.dt/6) * (k1 + 2*k2 + 2*k3 + k4);

            obj.get_chaser_pos();
            obj.get_chaser_euler();
        end
        
        % Dynamics (Internal Calculation)
        function dstate = dynamics(obj, x, u_ctrl, u_dist)
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
            
            w_t = obj.Target.stateECI(10:12);
            dw_t = zeros([3, 1]); % Temporal setting (placeholder)
            dv_t = obj.Target.gravitational_force();
            
            % Kinematics and Dynamics Setup
            R_tc = obj.get_Rt_c(s_curr);
            Rw_t = R_tc * w_t;
            w_c = w_curr + Rw_t; % Chaser relative angular velocity
            
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
            C2 = -Omega_wc;
            f_g = obj.gravitational_force();
            D2 = f_g - (R_tc * dv_t);
            dv = (C2 * v_curr) + D2 + ( f_ctrl + f_d ) / obj.m_c;
            
            dstate = [dsigma; domega; drho; dv];
        end
        
        function get_chaser_euler(obj)
            % GET_CHASER_EULER Computes Chaser's Euler angles from Target's Euler and Relative MRP
            %
            % Output:
            %   att_c: [roll, pitch, yaw] (rad) of Chaser w.r.t ECI
            
            phi = obj.Target.stateECI(7);
            theta = obj.Target.stateECI(8);
            psi = obj.Target.stateECI(9);
            
            R_i_t = angle2dcm(psi, theta, phi, 'ZYX');
            
            s = obj.state(1:3);
            s_sq = s' * s;
            skew_s = obj.skew(s);
            denom = (1 + s_sq)^2;
            
            R_t_c = eye(3) - (4*(1-s_sq)/denom)*skew_s + (8/denom)*(skew_s*skew_s);
            
            R_i_c = R_t_c * R_i_t;
            
            roll_c  = atan2(R_i_c(2,3), R_i_c(3,3));
            if R_i_c(1,3) > 1 || R_i_c(1, 3) < -1
                error("Out of range");
            end
            pitch_c = -asin(R_i_c(1,3));
            yaw_c   = atan2(R_i_c(1,2), R_i_c(1,1));
            
            obj.att_c = [roll_c; pitch_c; yaw_c];
        end

        function get_chaser_pos(obj)
            % GET_CHASER_POS Computes Chaser's position represented in ECI frame
            %
            % Output:
            %   pos_c: [x, y, z] (m) of Chaser w.r.t ECI
            
            obj.rho_c = obj.state(7:9) + obj.Target.stateECI(1:3);
        end

        function get_chaser_omg(obj)
            R_tc = obj.get_Rt_c(obj.state(1:3));
            w_t = obj.Target.stateECI(10:12);
            w_curr = obj.state(4:6);
            Rw_t = R_tc * w_t;
            obj.omg_c = w_curr + Rw_t; % Chaser relative angular velocity
        end
        
        function f_g = gravitational_force(obj)
            r_t = obj.Target.stateECI(1:3);
            R_tc = obj.get_Rt_c(obj.state(1:3));
            r_c = obj.state(7:9) + R_tc * r_t;
            r_c_norm = norm(r_c);
            f_g = -(obj.MU / r_c_norm^3) * r_c;
        end

        function c1 = get_C1(obj)
            w = obj.state(4:6);
            R_tc = obj.get_Rt_c(obj.state(1:3));
            w_t = obj.Target.stateECI(10:12);
            Rw_t = R_tc * w_t;
            w_c = w + Rw_t; % Chaser relative angular velocity

            S_Rwt = obj.skew(Rw_t);

            c1 = -obj.J_c*S_Rwt - S_Rwt*obj.J_c + obj.skew(obj.J_c*w_c);
        end

        function d1 = get_D1(obj)
            s = obj.state(1:3);
            R_tc = obj.get_Rt_c(s);
            w_t = obj.Target.stateECI(10:12);
            dw_t = zeros([3, 1]); % Temporal setting (placeholder)

            Rw_t = R_tc * w_t;

            S_Rwt = obj.skew(Rw_t);

            d1 = -S_Rwt*(obj.J_c*Rw_t) - obj.J_c*(R_tc*dw_t);
        end

        % Math Helpers
        function S = skew(~, x)
            S = [0,    -x(3),  x(2);
                 x(3),  0,    -x(1);
                -x(2),  x(1),  0];
        end
        
        function G = get_G_matrix(obj, s)
            G = 0.25 * ((1 - s' * s)*eye(3) + 2*obj.skew(s) + 2*(s*s'));
        end
        
        function R = get_Rt_c(obj, s)
            s_sq = s' * s;
            Omega_s = obj.skew(s);
            denom = (1 + s_sq)^2;
            R = eye(3) - (4*(1-s_sq)/denom)*Omega_s + (8/denom)*(Omega_s*Omega_s);
        end
    end
end