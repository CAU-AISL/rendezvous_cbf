classdef Logger < handle
    properties
        log
        index
    end

    methods
        function obj = Logger(logger_size)
            obj.log = zeros([logger_size.n, logger_size.m]);
            obj.index = 1;
        end

        function log_data(obj, data)
            obj.log(:, obj.index) = data;
            obj.index = obj.index + 1;
        end
    end
end