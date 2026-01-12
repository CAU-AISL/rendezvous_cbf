close all; clear all; clc;

addpath("utils\");

target_init_coe = struct('a', 7702455,...
                         'e', 0.12,...
                         'i', deg2rad(30),...
                         'Omega', deg2rad(0),...
                         'omega', deg2rad(0),...
                         'f0', deg2rad(0));

dt = 0.1;
time_span = 0:0.1:4000;

TargetSatellite = SatelliteDynamics(target_init_coe, dt);

target_state_log = zeros([length(TargetSatellite.stateECI), length(time_span)]);

for i = 1:length(time_span)
    
    TargetSatellite.step();
    target_state_log(:, i) = TargetSatellite.stateECI; 
end

run("plot_sim.m");
