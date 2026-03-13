classdef Logger < handle
    properties
        log
        index
    end

    methods
        function obj = Logger(n, m)
            obj.log = zeros([n, m]);
            obj.index = 1;
        end

        function log_data(obj, data)
            obj.log(:, obj.index) = data;
            obj.index = obj.index + 1;
        end
    end
end