classdef SatelliteDynamics < handle
    properties
        a       % Semi-major axis (COE)
        e       % Eccentricity (COE)
        i       % Inclination (COE)
        Omega   % Longitude of the ascending node (COE)
        omega   % Arguement of perigee (COE)
        f0      % Initial True anomaly (COE)

        MU      % Gravitational constant

        stateECI % State vector represented in ECI frame [position; velocity; attitude; angular rate] (12x1)
        stateCOE % State vector representing classical orbital elements (COEs)

        time    % Current simulation time [s]
        dt      % Time step
    end

    methods
        function obj = SatelliteDynamics(Cfg)
            init_coe = Cfg.target_init_state;
            obj.a       = init_coe.a;
            obj.e       = init_coe.e;
            obj.i       = init_coe.i;
            obj.Omega   = init_coe.Omega;
            obj.omega   = init_coe.omega;
            obj.f0      = init_coe.f0;

            obj.MU      = 3.986004e14;

            obj.stateECI = zeros(12, 1); % Initialize state vector in ECI frame
            obj.COE2ECI();
            obj.stateCOE = [obj.a;...
                obj.e;...
                obj.i;...
                obj.Omega;...
                obj.omega;...
                obj.f0];

            obj.time = 0;
            obj.dt = Cfg.dt;
        end

        function step(obj)
            x_curr = obj.stateECI;

            % RK4 Integration
            k1 = obj.dynamics(x_curr);
            k2 = obj.dynamics(x_curr + obj.dt/2*k1);
            k3 = obj.dynamics(x_curr + obj.dt/2*k2);
            k4 = obj.dynamics(x_curr + obj.dt*k3);

            % Update the current time and state vector
            obj.time = obj.time + obj.dt;
            obj.stateECI = x_curr + (obj.dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        end

        function dstate = dynamics(obj, x)
            pos = x(1:3);
            vel = x(4:6);
            euler = x(7:9); % [phi; theta; psi] (Roll, Pitch, Yaw)
            omg = x(10:12);
            
            r = norm(pos);
            
            acc = -obj.MU * pos / r^3;
            
            % Z-Y-X Sequence (Standard Aerospace)
            phi = euler(1); 
            theta = euler(2);

            c_phi = cos(phi);
            s_phi = sin(phi);
            c_th  = cos(theta);
            s_th  = sin(theta);

            if abs(c_th) < 1e-4
                c_th = sign(c_th) * 1e-4;
                if c_th == 0, c_th = 1e-4; end
                warning("Singularity approaching at theta = %.2f deg", rad2deg(theta));
            end

            t_th = s_th / c_th;

            H = [1,  s_phi * t_th,   c_phi * t_th;
                0,  c_phi,         -s_phi;
                0,  s_phi / c_th,   c_phi / c_th];
             
            d_euler = H * omg;
            
            % Assuming Target is a rigid body with constant angular velocity 
            % (Torque-free motion approximation for simple target)
            d_omg = zeros(3, 1); 
            % If you have inertia J, use: d_omg = J \ (-cross(omg, J*omg));
            
            dstate = [vel; acc; d_euler; d_omg];
        end

        function COE2ECI(obj)
            % h: Angular momentum
            h = sqrt(obj.MU*obj.a*(1 - obj.e^2));
            rotm_P2E = obj.rotm_PFC2ECI();

            % Position in perifocal frame
            p = (h^2 / obj.MU) * (1 / (1 + obj.e * cos(obj.f0))) * [cos(obj.f0); sin(obj.f0); 0];
            obj.stateECI(1:3) = rotm_P2E*p;

            % Velocity in perifocal frame
            v = (obj.MU / h) * [-sin(obj.f0); obj.e + cos(obj.f0); 0];
            obj.stateECI(4:6) = rotm_P2E*v;
        end

        function R = rotm_PFC2ECI(obj)
            % Define rotation matrix from ECI to perifocal frame
            cosOmg = cos(obj.Omega);
            sinOmg = sin(obj.Omega);
            cosomg = cos(obj.omega);
            sinomg = sin(obj.omega);
            cosInc = cos(obj.i);
            sinInc = sin(obj.i);

            RzOmg = [ cosOmg, sinOmg, 0;...
                -sinOmg, cosOmg, 0;...
                0,      0, 1];
            RxInc = [1,       0,      0;...
                0,  cosInc, sinInc;...
                0, -sinInc, cosInc];
            Rzomg = [cosomg, sinomg, 0;...
                -sinomg, cosomg, 0;...
                0,      0, 1];

            R = (Rzomg*RxInc*RzOmg)';
        end

        function f_tb = gravitational_force(obj)
            % Compute gravitational force acting on the target satellite in target frame
            r_t = obj.stateECI(1:3);
            r = norm(r_t);
            f_t = -obj.MU * r_t / r^3;
            att = obj.stateECI(7:9);
            f_tb = angle2dcm(att(3), att(2), att(1), 'ZYX') * f_t;
        end
    end
end