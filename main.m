close all; clear all; clc;
rng(0);
addpath(genpath(pwd));

simCfg = SimCfg();
simLogger = LoggerCfg(simCfg.sim_len);

TargetSatellite = SatelliteDynamics(simCfg);

ChaserSatellite = RelativeDynamics(simCfg, TargetSatellite);
% ChaserSatellite = StochasticRelative(simCfg, TargetSatellite);

% Controller = ClfQp(simCfg, ChaserSatellite);
% Controller = HOCBF(simCfg, ChaserSatellite);
Controller = CCBF(simCfg, ChaserSatellite);
% Controller = SCCBF(simCfg, ChaserSatellite);

u_dist.tau_d = zeros([3, 1]);
u_dist.f_d = zeros([3, 1]);

for i = simCfg.sim_time
    u_ctrl = Controller.command();
    
    TargetSatellite.step();
    ChaserSatellite.step(u_ctrl, u_dist);

    simLogger.loggerTarget.state.log_data(TargetSatellite.stateECI);
    simLogger.loggerRelative.state.log_data(ChaserSatellite.state);
    simLogger.loggerChaser.pos.log_data(ChaserSatellite.rho_c);
    simLogger.loggerChaser.att.log_data(ChaserSatellite.att_c);
    simLogger.loggerControl.vel_d.log_data(Controller.ref_vel);
    simLogger.loggerControl.force.log_data(u_ctrl.f);
    simLogger.loggerControl.omg_d.log_data(Controller.ref_omg);
    simLogger.loggerControl.moment.log_data(u_ctrl.tau);
    simLogger.loggerLyapunov.V.log_data([Controller.rhoV; Controller.velV; Controller.sigV; Controller.omgV]);
    simLogger.loggerBarrier.h.log_data(Controller.barrier_value());
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
saveFilename = fullfile(pwd, sprintf('simLogger_%s.mat', timestamp));
save(saveFilename, 'simLogger');

run("plot_sim.m");
