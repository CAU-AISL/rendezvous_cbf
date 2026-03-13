close all;
% Ready for plotting the Earth

plot_earth = false;
if plot_earth
    [x, y, z] = sphere(50);
    Re = 6371;
    x = x*Re;
    y = y*Re;
    z = z*Re;
    load('topo.mat', 'topo'); 
    orbit_width = 1.5;
    
    earthPlot = figure();
    earthPlot.Theme = 'light';
    target_orbit_m2km = TargetLogger.state.log(1:3, :).*(10^-3);
    chaser_orbit_m2km = ChaserLogger.pos.log(:, :).*(10^-3);
    hold on; grid on;
    plot3(target_orbit_m2km(1,:), target_orbit_m2km(2,:), target_orbit_m2km(3,:), 'LineWidth', orbit_width);
    plot3(chaser_orbit_m2km(1,:), chaser_orbit_m2km(2,:), chaser_orbit_m2km(3,:), 'LineWidth', orbit_width);
    s = surface(x, y, z);
    s.FaceColor = 'texturemap';    % Activate texture mapping
    s.CData = topo;                % Assign terrain data
    s.EdgeColor = 'none';          % Remove grid lines
    s.FaceAlpha = 0.3;
    axis equal;
    view(3);
    colormap([zeros(64,1), zeros(64,1), linspace(0.5,1,64)';...
                  linspace(0.5,0.8,64)', linspace(1,0.7,64)', zeros(64,1)]);
    xlabel("x-axis (km)");
    ylabel("y-axis (km)");
    zlabel("z-axis (km)");
end
%%
fig_size = [800, 600];
line_width = 1.5;

posPlot = figure();
posPlot.Theme = 'light';
posPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
pos = plot(simCfg.sim_time, simLogger.loggerRelative.state.log(7,:), 'LineWidth', line_width);
rfp = plot(simCfg.sim_time, ones(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\rho_x (m)', 'Interpreter', 'tex');
legend([pos, rfp], {'Response', 'Reference'});
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerRelative.state.log(8,:), 'LineWidth', line_width);
plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\rho_y (m)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerRelative.state.log(9,:), 'LineWidth', line_width);
plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\rho_z (m)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/pos_plot.png');

velPlot = figure();
velPlot.Theme = 'light';
velPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
vel = plot(simCfg.sim_time, simLogger.loggerRelative.state.log(10,:), 'LineWidth', line_width);
% rfv = plot(simCfg.sim_time, simLogger.loggerControl.vel_d.log(1,:), 'r--', 'LineWidth', line_width);
ylabel('V_x (m/s)', 'Interpreter', 'tex');
% legend([vel, rfv], {'Response', 'Reference'}, 'Location', 'southeast');
legend([vel], {'Response'}, 'Location', 'northeast');
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerRelative.state.log(11,:), 'LineWidth', line_width);
% plot(simCfg.sim_time, simLogger.loggerControl.vel_d.log(2,:), 'r--', 'LineWidth', line_width);
ylabel('V_y (m/s)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerRelative.state.log(12,:), 'LineWidth', line_width);
% plot(simCfg.sim_time, simLogger.loggerControl.vel_d.log(3,:), 'r--', 'LineWidth', line_width);
ylabel('V_z (m/s)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/vel_plot.png');

attPlot = figure();
attPlot.Theme = 'light';
attPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
att = plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(1,:)), 'LineWidth', line_width);
rfa = plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\sigma_x (deg)', 'Interpreter', 'tex');
legend([att, rfa], {'Response', 'Reference'});
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(2,:)), 'LineWidth', line_width);
plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\sigma_y (deg)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(3,:)), 'LineWidth', line_width);
plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylabel('\sigma_z (deg)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/att_plot.png');

omgPlot = figure();
omgPlot.Theme = 'light';
omgPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
omg = plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(4,:)), 'LineWidth', line_width);
% rfo = plot(simCfg.sim_time, rad2deg(simLogger.loggerControl.omg_d.log(1,:)), 'r--', 'LineWidth', line_width);
% legend([omg, rfo], {'Response', 'Reference'});
legend([omg], {'Response'}, 'Location', 'northeast');
ylabel('\omega_x (deg/s)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(5,:)), 'LineWidth', line_width);
% plot(simCfg.sim_time, rad2deg(simLogger.loggerControl.omg_d.log(2,:)), 'r--', 'LineWidth', line_width);
ylabel('\omega_y (deg/s)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, rad2deg(simLogger.loggerRelative.state.log(6,:)), 'LineWidth', line_width);
% plot(simCfg.sim_time, rad2deg(simLogger.loggerControl.omg_d.log(3,:)), 'r--', 'LineWidth', line_width);
ylabel('\omega_z (deg/s)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/omg_plot.png');

forcePlot = figure();
forcePlot.Theme = 'light';
forcePlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.force.log(1,:), 'LineWidth', line_width);
yline([-20, 20], 'r--');
ylim([-22, 22]);
ylabel('F_x (N)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.force.log(2,:), 'LineWidth', line_width);
yline([-20, 20], 'r--');
ylim([-22, 22]);
ylabel('F_y (N)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.force.log(3,:), 'LineWidth', line_width);
yline([-20, 20], 'r--');
ylim([-22, 22]);
ylabel('F_z (N)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/force_plot.png');

momentPlot = figure();
momentPlot.Theme = 'light';
momentPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.moment.log(1,:), 'LineWidth', line_width);
yline([-5, 5], 'r--');
ylim([-6, 6]);
ylabel('M_x (N)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.moment.log(2,:), 'LineWidth', line_width);
yline([-5, 5], 'r--');
ylim([-6, 6]);
ylabel('M_y (N)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerControl.moment.log(3,:), 'LineWidth', line_width);
yline([-5, 5], 'r--');
ylim([-6, 6]);
ylabel('M_z (N)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/moment_plot.png');

lyapunovPlot = figure();
lyapunovPlot.Theme = 'light';
lyapunovPlot.Position(3:4) = fig_size;
subplot(4, 1, 1);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerLyapunov.V.log(1,:), 'LineWidth', line_width);
ylabel('V_{\rho}', 'Interpreter', 'tex');
subplot(4, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerLyapunov.V.log(2,:), 'LineWidth', line_width);
ylabel('V_v', 'Interpreter', 'tex');
subplot(4, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerLyapunov.V.log(3,:), 'LineWidth', line_width);
ylabel('V_{\sigma}', 'Interpreter', 'tex');
subplot(4, 1, 4);
hold on; grid on;
plot(simCfg.sim_time, simLogger.loggerLyapunov.V.log(4,:), 'LineWidth', line_width);
ylabel('V_{\omega}', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/lyapunov_plot.png');

barrierPlot = figure();
barrierPlot.Theme = 'light';
barrierPlot.Position(3:4) = fig_size;
hold on; grid on;
hPlot1 = plot(simCfg.sim_time, simLogger.loggerBarrier.h.log(1,:), 'LineWidth', line_width);
hPlot2 = plot(simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r-', 'LineWidth', line_width);
uistack(hPlot1, 'top');
ylabel('h', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/barrier_plot.png');

rt = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    sigma = simLogger.loggerRelative.state.log(1:3, i);
    rho = simLogger.loggerRelative.state.log(7:9, i);
    rt(:, i) = (ChaserSatellite.get_R_tc(sigma))'*rho;
end

offset = 1;
x_min = 1;
x_max = max(rt(1,:)) + offset;
x = linspace(x_min, x_max, 100);
theta = linspace(0, 2*pi, 50);
[X, Theta] = meshgrid(x, theta);
R = sqrt(controlCfg.a_h*(X - controlCfg.delta_h).^3);
Y = R .* cos(Theta);
Z = R .* sin(Theta);

rtPlot = figure();
rtPlot.Theme = 'light';
rtPlot.Position(3:4) = fig_size;
hold on; grid on;
plot3(rt(1, :), rt(2, :), rt(3, :), 'LineWidth', line_width);
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'r');

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

view(3);
xlabel('x_t [m]'); ylabel('y_t [m]'); zlabel('z_t [m]');
xlim([-5, 30]); ylim([-20, 20]); zlim([-20, 20]);
saveas(gcf, 'assets/rt_plot.png');

%%
sim_flag = false;
% sim_flag = true;
if sim_flag
    epochYear = 2026;
    epochMonth = 1;
    epochDay = 9;
    epochHour = 0;
    epochMinute = 0;
    epochSecond = 0;
    epoch = datetime(epochYear, epochMonth, epochDay, epochHour, epochMinute, epochSecond);
    pos_name = '<Xicrf>';
    att_name = '<qicrf2b>';
    
    startTime = epoch;
    stopTime = startTime + seconds(simCfg.sim_time(end));
    sc = satelliteScenario(startTime,stopTime,simCfg.dt);
    targetPos = timeseries(simLogger.loggerTarget.state.log(1:3, :)', simCfg.sim_time);
    targetPos.Name = pos_name;
    targetAtt = timeseries(eul2quat(simLogger.loggerTarget.state.log(7:9, :)'), simCfg.sim_time);
    targetAtt.Name = att_name;
    target_sat = satellite(sc, targetPos, Name='Target_Sat');
    pointAt(target_sat, targetAtt);
    
    chaserPos = timeseries(simLogger.loggerChaser.pos.log(1:3, :)', simCfg.sim_time);
    chaserPos.Name = pos_name;
    chaserAtt = timeseries(eul2quat(simLogger.loggerChaser.att.log(1:3, :)'), simCfg.sim_time);
    chaserAtt.Name = att_name;
    chaser_sat = satellite(sc, chaserPos, Name='Chaser_Sat');
    pointAt(chaser_sat, chaserAtt);
    
    ScenarioViewer = satelliteScenarioViewer(sc, CameraReferenceFrame='Inertial');
    ScenarioViewer.Position(3:4) = [1600, 900];
    ScenarioViewer.PlaybackSpeedMultiplier = 3;
    target_sat.Visual3DModel = 'SmallSat.glb';
    target_sat.Visual3DModelScale = 0.8;
    chaser_sat.Visual3DModel = 'SmallSat.glb';
    chaser_sat.Visual3DModelScale = 0.8;
    hide(target_sat.Orbit);
    hide(chaser_sat.Orbit);
    play(sc);
    camtarget(ScenarioViewer, target_sat);
end

