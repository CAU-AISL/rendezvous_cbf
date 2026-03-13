close all; clear all; clc;

addpath(genpath(pwd));

fig_size = [600, 450];
line_width = 1.5;

load('assets\nominal_clfqp\simLogger_20260312_184550.mat');
nominal = struct('sig', simLogger.loggerRelative.state.log(1:3, :),...
                 'rho', simLogger.loggerRelative.state.log(7:9, :));


load('assets\ccbf\simLogger_20260313_125027.mat');
ccbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
              'force', simLogger.loggerControl.force.log,...
              'sig', simLogger.loggerRelative.state.log(1:3, :),...
              'rho', simLogger.loggerRelative.state.log(7:9, :));

load('assets\hocbf\simLogger_20260313_125227.mat');
hocbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
               'force', simLogger.loggerControl.force.log,...
               'sig', simLogger.loggerRelative.state.log(1:3, :),...
               'rho', simLogger.loggerRelative.state.log(7:9, :));

simCfg = SimCfg();
controlCfg = ControlCfg();

lyapunovFig = figure();
lyapunovFig.Theme = 'light';
lyapunovFig.Position(3:4) = fig_size;
subplot(2, 1, 1);
hold on; grid on;
CcbfLyp = plot(simCfg.sim_time, ccbf.lyapunov(1,:), 'LineWidth', line_width);
HocbfLyp = plot(simCfg.sim_time, hocbf.lyapunov(1,:), 'LineWidth', line_width);
ylabel('V_\rho'); legend([CcbfLyp, HocbfLyp], {'Cascaded CBF', 'HOCBF'}, 'Location', 'northeast');

subplot(2, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, ccbf.lyapunov(2,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.lyapunov(2,:), 'LineWidth', line_width);
ylabel('V_v');

inputFig = figure();
inputFig.Theme = 'light';
inputFig.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(1,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(1,:), 'LineWidth', line_width);
ylabel('F_x');
legend({'Cascaded CBF', 'HOCBF'}, 'Location', 'northeast');

subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(2,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(2,:), 'LineWidth', line_width);
ylabel('F_y');

subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(3,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(3,:), 'LineWidth', line_width);
ylabel('F_z');

nominal_traj = NaN([3, simCfg.sim_len]);
ccbf_traj = NaN([3, simCfg.sim_len]);
hocbf_traj = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    nominal_sigma = nominal.sig(:, i);
    nominal_rho = nominal.rho(:, i);
    nominal_traj(:, i) = (get_R_tc(nominal_sigma))'*nominal_rho;

    ccbf_sigma = ccbf.sig(:, i);
    ccbf_rho = ccbf.rho(:, i);
    ccbf_traj(:, i) = (get_R_tc(ccbf_sigma))'*ccbf_rho;

    hocbf_sigma = hocbf.sig(:, i);
    hocbf_rho = hocbf.rho(:, i);
    hocbf_traj(:, i) = (get_R_tc(hocbf_sigma))'*hocbf_rho;
end

rtPlot = figure();
rtPlot.Theme = 'light';
rtPlot.Position(3:4) = [947, 506];
hold on; grid on;
plot3(ccbf_traj(1, :), ccbf_traj(2, :), ccbf_traj(3, :), 'LineWidth', line_width);
plot3(hocbf_traj(1, :), hocbf_traj(2, :), hocbf_traj(3, :), 'LineWidth', line_width);
plot3(nominal_traj(1, :), nominal_traj(2, :), nominal_traj(3, :), 'LineWidth', line_width);
% Surface plotting removed as requested

meshData = readSurfaceMesh('SmallSat.glb');

V = double(meshData.Vertices);
F = double(meshData.Faces);
V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';

scale_factor = 0.8; 
V = V * scale_factor;

trisurf(F, V(:,1), V(:,2), V(:,3), ...
      'FaceColor', [0.7 0.7 0.7], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');

camlight('headlight');
legend({'Cascaded CBF', 'HOCBF', 'Nominal Controller'}, 'Location', 'best');

% view(3);
view([-26.5, 25.5])
xlabel('x_t [m]'); ylabel('y_t [m]'); zlabel('z_t [m]');
% xlim([-5, 30]); ylim([-20, 20]); zlim([-20, 20]);
axis equal;

%%

function R = get_R_tc(mrp)
    mrp_sq = mrp' * mrp;
    skew = [      0, -mrp(3),  mrp(2);...
             mrp(3),       0, -mrp(1);...
            -mrp(2),  mrp(1),       0];
    denom = (1 + mrp_sq)^2;
    R = eye(3) - (4*(1-mrp_sq)/denom)*skew + (8/denom)*(skew*skew);
end
