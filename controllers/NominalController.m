classdef NominalController < handle
    % Nominal Controller Base Class
    properties
        settingCtrlRho
        settingCtrlVel
        settingCtrlSig
        settingCtrlOmg

        RelativeChaser  RelativeDynamics
        nominalRho  NominalRho
        nominalVel  NominalVel
        nominalSig  NominalSigma
        nominalOmg  NominalOmega

        nominalLogRho
        nominalLogVel
        nominalLogSig
        nominalLogOmg
    end

    methods
        function obj = NominalController(relativeDynamics)
            obj.settingCtrlRho = struct('gamma', 0.5, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e6, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.settingCtrlVel = struct('gamma', 0.5, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e12, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.settingCtrlSig = struct('gamma', 0.2, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e12, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.settingCtrlOmg = struct('gamma', 0.2, ...
                                        'alpha', 1.0, ...
                                        'p_relax', 1e12, ...
                                        'qp_options', optimoptions('quadprog', 'Display', 'off'));
            obj.RelativeChaser = relativeDynamics;
            obj.nominalRho = NominalRho(obj.settingCtrlRho, obj.RelativeChaser);
            obj.nominalVel = NominalVel(obj.settingCtrlVel, obj.RelativeChaser);
            obj.nominalSig = NominalSigma(obj.settingCtrlSig, obj.RelativeChaser);
            obj.nominalOmg = NominalOmega(obj.settingCtrlOmg, obj.RelativeChaser);

            obj.nominalLogRho = struct('Vd', zeros([3, 1]),...
                                        'slack', 0,...
                                        'feas', 0);
            obj.nominalLogVel = struct('Force', zeros([3, 1]),...
                                        'slack', 0,...
                                        'feas', 0);
            obj.nominalLogSig = struct('Omgd', zeros([3, 1]),...
                                        'slack', 0,...
                                        'feas', 0);
            obj.nominalLogOmg = struct('Moment', zeros([3, 1]),...
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
            % u = Force;
            u.f = Force;

            [omg_d, slack_sig, feas_sig] = obj.nominalSig.command();
            obj.nominalLogSig.Omgd = omg_d;
            obj.nominalLogSig.slack = slack_sig;
            obj.nominalLogSig.feas = feas_sig;

            [Moment, slack_omg, feas_omg] = obj.nominalOmg.command(omg_d);
            obj.nominalLogOmg.Moment = Moment;
            obj.nominalLogOmg.slack = slack_omg;
            obj.nominalLogOmg.feas = feas_omg;
            % u = Moment;
            u.tau = Moment;
        end
    end
end

