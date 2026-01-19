classdef NominalController < handle
    % Nominal Controller Base Class
    properties
        settingCtrlRho
        settingCtrlVel
        RelativeChaser  RelativeDynamics
        nominalRho  NominalRho
        nominalVel  NominalVel
        nominalLogRho
        nominalLogVel
    end

    methods
        function obj = NominalController(relativeDynamics)
            obj.settingCtrlRho = struct('gamma', 5.0, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e6, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.settingCtrlVel = struct('gamma', 0.5, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e12, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.RelativeChaser = relativeDynamics;
            obj.nominalRho = NominalRho(obj.settingCtrlRho, obj.RelativeChaser);
            obj.nominalVel = NominalVel(obj.settingCtrlVel, obj.RelativeChaser);
            obj.nominalLogRho = struct('Vd', zeros([3, 1]),...
                                        'slack', 0,...
                                        'feas', 0);
            obj.nominalLogVel = struct('Force', zeros([3, 1]),...
                                        'slack', 0,...
                                        'feas', 0);
        end

        function u = compute_control(obj)
            % Compute desired rho_dot
            [Vd, slack_rho, feas_rho] = obj.nominalRho.command();
            obj.nominalLogRho.Vd = Vd;
            obj.nominalLogRho.slack = slack_rho;
            obj.nominalLogRho.feas = feas_rho;

            % Update desired state for velocity controller
            [Force, slack_vel, feas_vel] = obj.nominalVel.command(obj.nominalLogRho.Vd);
            obj.nominalLogVel.Force = Force;
            obj.nominalLogVel.slack = slack_vel;
            obj.nominalLogVel.feas = feas_vel;
            u = Force;
        end
    end
end

