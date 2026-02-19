close all; clear; clc;

addpath(genpath(pwd));

simCfg = SimCfg();
simLogger = LoggerCfg(simCfg.sim_len);
controlCfg = ControlCfg();

TargetSatellite = SatelliteDynamics(simCfg.target_init_state, simCfg.dt);
ChaserSatellite = RelativeDynamics(simCfg.chaser_init_state, simCfg.dt, TargetSatellite);
Controller = ClfQp(controlCfg, ChaserSatellite);

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
end

run("plot_sim.m");
