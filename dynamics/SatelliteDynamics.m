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
        function obj = SatelliteDynamics(init_coe, dt)
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
            obj.dt = dt;
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
            att = x(7:9);
            omg = x(10:12);

            % Compute the gravitational acceleration
            r = norm(pos);
            acc = -obj.MU * pos / r^3;

            % Update the state derivative
            dstate = zeros(12, 1);
            dstate(1:3) = vel;  % Position derivative
            dstate(4:6) = acc;  % Velocity derivative
            dstate(7:9) = zeros(3, 1);  % Attitude dynamics (placeholder)
            dstate(10:12) = zeros(3, 1); % Angular rate dynamics (placeholder)
        end

        function COE2ECI(obj)
            % h: Angular momentum
            h = sqrt(obj.MU*obj.a*(1 - obj.e^2));
            rotm_E2P = obj.rotm_ECI2PFC();

            % Position in perifocal frame
            p = (h^2 / obj.MU) * (1 / (1 + obj.e * cos(obj.f0))) * [cos(obj.f0); sin(obj.f0); 0];
            obj.stateECI(1:3) = rotm_E2P*p;
            
            % Velocity in perifocal frame
            v = (obj.MU / h) * [-sin(obj.f0); obj.e + cos(obj.f0); 0];
            obj.stateECI(4:6) = rotm_E2P*v;
        end

        function R = rotm_ECI2PFC(obj)
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
    end
end