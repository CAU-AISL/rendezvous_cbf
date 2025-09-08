close all; clear; clc;

addpath("dynamics\");
addpath("controllers\");
addpath("utils\");

target_init_coe = struct('a', 7702455,...
                         'e', 0.12,...
                         'i', deg2rad(30),...
                         'Omega', deg2rad(0),...
                         'omega', deg2rad(0),...
                         'f0', deg2rad(0));

relative_init_state = struct('sigma', [-0.1; -0.2; 0.1],...
                             'omega', [0.1; -0.2; 0.4],...
                             'rho', [30; 10; -20],...
                             'vel', [0.5; -0.5; 0.1]);

dt = 0.01;
time_span = 0:dt:100;
sim_length = length(time_span);
TargetSatellite = SatelliteDynamics(target_init_coe, dt);
state_log_size = struct('n', 12,...
                       'm', sim_length);

TargetLogger = struct('state', Logger(state_log_size));

ChaserSatellite = RelativeDynamics(relative_init_state, dt, TargetSatellite);
RelativeLogger = struct('state', Logger(state_log_size));

chaser_state.n = 3;
chaser_state.m = sim_length;
ChaserLogger = struct('pos', Logger(chaser_state),...
                      'att', Logger(chaser_state));

target_state_log = zeros([length(TargetSatellite.stateECI), length(time_span)]);

for i = 1:length(time_span)
    u_ctrl.tau = zeros([3, 1]);
    u_ctrl.f = zeros([3, 1]);
    u_dist.tau_d = zeros([3, 1]);
    u_dist.f_d = zeros([3, 1]);
    TargetSatellite.step();
    ChaserSatellite.step(u_ctrl, u_dist);

    TargetLogger.state.log_data(TargetSatellite.stateECI);
    RelativeLogger.state.log_data(ChaserSatellite.state);
    ChaserLogger.pos.log_data(ChaserSatellite.rho_c);
    ChaserLogger.att.log_data(ChaserSatellite.att_c);
end

run("plot_sim.m");
