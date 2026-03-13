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
            obj.T   = 200;
            obj.sim_time = 0 : obj.dt : obj.T;
            obj.sim_len = length(obj.sim_time);
            obj.target_init_state = struct('a', 7702455,...
                                           'e', 0.12,...
                                           'i', deg2rad(30),...
                                           'Omega', deg2rad(0),...
                                           'omega', deg2rad(0),...
                                           'f0', deg2rad(0));
            obj.chaser_init_state =  struct('sigma', [-0.1; 0.12; 0.1],...
                                            'omega', [0.05; -0.03; 0.07],...
                                            'rho', [23.6; -8.3; 19.2],...
                                            'vel', [-0.2; -0.3; -0.1]);
        end
    end
end