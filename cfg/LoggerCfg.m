classdef LoggerCfg
    properties
        loggerTarget
        loggerRelative
        loggerChaser
        loggerControl
        loggerLyapunov
        loggerBarrier
    end

    methods
        function obj = LoggerCfg(sim_len)
            obj.loggerTarget = struct('state', Logger(12, sim_len));
            obj.loggerRelative = struct('state', Logger(12, sim_len));
            obj.loggerChaser = struct('pos', Logger(3, sim_len),...
                                      'att', Logger(3, sim_len));
            obj.loggerControl =  struct('force', Logger(3, sim_len),...
                                        'vel_d', Logger(3, sim_len),...
                                        'omg_d', Logger(3, sim_len),...
                                        'moment', Logger(3, sim_len));
            obj.loggerLyapunov = struct('V', Logger(4, sim_len));
            obj.loggerBarrier = struct('h', Logger(1, sim_len));
        end
    end
end