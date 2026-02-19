classdef SimCfg < handle
    properties
        dt
        T
        sim_time
        sim_len
        target_init_state
        chaser_init_state
    end
    methods
        function obj = SimCfg()
            obj.dt  = 0.01;
            obj.T   = 100;
            obj.sim_time = 0 : obj.dt : obj.T;
            obj.sim_len = length(obj.sim_time);
            obj.target_init_state = struct('a', 7702455,...
                                           'e', 0.12,...
                                           'i', deg2rad(30),...
                                           'Omega', deg2rad(0),...
                                           'omega', deg2rad(0),...
                                           'f0', deg2rad(0));
            obj.chaser_init_state =  struct('sigma', [-0.1; -0.2; 0.1],...
                                            'omega', [0.1; -0.2; 0.4],...
                                            'rho', [30; 10; -20],...
                                            'vel', [0.5; -0.5; 0.1]);
        end
    end
end