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
fig_size = [800, 450];
line_width = 1.5;

posPlot = figure();
posPlot.Theme = 'light';
posPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(7,:), 'LineWidth', line_width);
ylabel('\rho_x (m)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(8,:), 'LineWidth', line_width);
ylabel('\rho_y (m)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(9,:), 'LineWidth', line_width);
ylabel('\rho_z (m)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/pos_plot.png');

velPlot = figure();
velPlot.Theme = 'light';
velPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(10,:), 'LineWidth', line_width);
ylabel('V_x (m/s)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(11,:), 'LineWidth', line_width);
ylabel('V_y (m/s)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(time_span, RelativeLogger.state.log(12,:), 'LineWidth', line_width);
ylabel('V_z (m/s)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/vel_plot.png');

attPlot = figure();
attPlot.Theme = 'light';
attPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(1,:)), 'LineWidth', line_width);
ylabel('\sigma_x (deg)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(2,:)), 'LineWidth', line_width);
ylabel('\sigma_y (deg)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(3,:)), 'LineWidth', line_width);
ylabel('\sigma_z (deg)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/att_plot.png');

omgPlot = figure();
omgPlot.Theme = 'light';
omgPlot.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(4,:)), 'LineWidth', line_width);
ylabel('\omega_x (deg/s)', 'Interpreter', 'tex');
subplot(3, 1, 2);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(5,:)), 'LineWidth', line_width);
ylabel('\omega_y (deg/s)', 'Interpreter', 'tex');
subplot(3, 1, 3);
hold on; grid on;
plot(time_span, rad2deg(RelativeLogger.state.log(6,:)), 'LineWidth', line_width);
ylabel('\omega_z (deg/s)', 'Interpreter', 'tex');
xlabel('Time (s)', 'Interpreter', 'tex');
saveas(gcf, 'assets/omg_plot.png');

%%
sim_flag = false;
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
    stopTime = startTime + seconds(time_span(end));
    sc = satelliteScenario(startTime,stopTime,dt);
    targetPos = timeseries(TargetLogger.state.log(1:3, :)', time_span);
    targetPos.Name = pos_name;
    targetAtt = timeseries(eul2quat(TargetLogger.state.log(7:9, :)'), time_span);
    targetAtt.Name = att_name;
    target_sat = satellite(sc, targetPos, Name='Target_Sat');
    pointAt(target_sat, targetAtt);
    
    chaserPos = timeseries(ChaserLogger.pos.log(1:3, :)', time_span);
    chaserPos.Name = pos_name;
    chaserAtt = timeseries(eul2quat(ChaserLogger.att.log(1:3, :)'), time_span);
    chaserAtt.Name = att_name;
    chaser_sat = satellite(sc, chaserPos, Name='Chaser_Sat');
    pointAt(chaser_sat, chaserAtt);
    
    ScenarioViewer = satelliteScenarioViewer(sc, CameraReferenceFrame='Inertial');
    target_sat.Visual3DModel = 'SmallSat.glb';
    target_sat.Visual3DModelScale = 0.8;
    chaser_sat.Visual3DModel = 'SmallSat.glb';
    chaser_sat.Visual3DModelScale = 0.8;
    % hide(target_sat.Orbit);
    play(sc);
    camtarget(ScenarioViewer, target_sat);
end

