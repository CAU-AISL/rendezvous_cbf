close all; clear; clc;
rng(0);
addpath(genpath(pwd));

N_mc = 500;
load("MonteCarloInitCondition.mat");

p = gcp('nocreate');
if isempty(p)
    parpool('local'); 
end

mc_loggers = cell(N_mc, 1);

tic;

parfor mc_idx = 1:N_mc

    simCfg = SimCfg();
    
    simCfg.chaser_init_state = chaser_init_state(mc_idx);
    
    TargetSatellite = SatelliteDynamics(simCfg);
    
    ChaserSatellite = RelativeDynamics(simCfg, TargetSatellite);
    
    Controller = HOCBF(simCfg, ChaserSatellite);
    % Controller = CCBF(simCfg, ChaserSatellite);
    
    simLogger = LoggerCfg(simCfg.sim_len);
    
    u_dist = struct();
    u_dist.tau_d = zeros(3, 1);
    u_dist.f_d = zeros(3, 1);
    
    
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
    
    mc_loggers{mc_idx} = simLogger;
end

time_elapsed = toc;
disp(['Monte-Carlo Simulation Conducted! Time elapsed: ', num2str(time_elapsed), 's']);

save('MC_Results_HOCBF.mat', 'mc_loggers', 'mc_init_states');
% save('MC_Results_CCBF.mat', 'mc_loggers', 'mc_init_states');

