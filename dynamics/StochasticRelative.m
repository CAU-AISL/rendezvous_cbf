classdef StochasticRelative < RelativeDynamics

    properties
        % sigma denotes standarad deviation of Normal Distribution
        sig_rho
        sig_vel
        sig_mrp
        sig_omg

        eps_rho
        eps_vel
        eps_mrp
        eps_omg
    end

    methods
        function obj = StochasticRelative(initial_state, dt, targetSatellite)
            obj@RelativeDynamics(initial_state, dt, targetSatellite);
            
            obj.sig_rho = 0.05;
            obj.sig_vel = 0.01;
            obj.sig_mrp = deg2rad(1);
            obj.sig_omg = deg2rad(0.5);

            % obj.sig_rho = 0.005;
            % obj.sig_vel = 0.01;
            % obj.sig_mrp = 0.002;
            % obj.sig_omg = 0.0005;

            obj.eps_rho = 0;
            obj.eps_vel = 0;
            obj.eps_mrp = 0;
            obj.eps_omg = 0;
        end

        function dstate = dynamics(obj, x, u_ctrl, u_dist)
            dstate_ = dynamics@RelativeDynamics(obj, x, u_ctrl, u_dist);
            dstate = dstate_ + [normrnd(obj.eps_mrp, obj.sig_mrp, 3, 1);...
                                normrnd(obj.eps_omg, obj.sig_omg, 3, 1);...
                                normrnd(obj.eps_rho, obj.sig_rho, 3, 1);...
                                normrnd(obj.eps_vel, obj.sig_vel, 3, 1)];
        end

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
            obj.R_tc = obj.get_R_tc(obj.state(1:3));
        end
    end
end