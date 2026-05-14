classdef LPF < handle
    properties
        omg_c   % Cut-off Frequency (rad/s)
        T       % Sampling time (s)
        alpha   % 1st-order LPF coefficient
        prev_u  % Previous output
    end
    methods
        function obj = LPF(omg_c, Cfg)
            obj.omg_c = omg_c;
            obj.T = Cfg.dt;
            obj.alpha = exp(-obj.omg_c * obj.T);
            obj.prev_u = [];
        end
        function filtered_u = forward(obj, u)
            if isempty(obj.prev_u)
                obj.prev_u = zeros(size(u));
            end
            filtered_u = (1 - obj.alpha) * u + obj.prev_u * obj.alpha;
            obj.prev_u = filtered_u;
        end
    end
end