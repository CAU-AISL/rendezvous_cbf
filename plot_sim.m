% Ready for plotting the Earth
[x, y, z] = sphere(50);
Re = 6371;
x = x*Re;
y = y*Re;
z = z*Re;
load('topo.mat', 'topo'); 
orbit_width = 1.5;

posPlot = figure();
posPlot.Theme = 'light';
orbit_m2km = target_state_log(1:3, :)./1000;
hold on; grid on;
% target = readSurfaceMesh("SmallSat.glb");
% sigma = [-0.1; -0.2; 0.1];
% rotate(target, "euler", rad2deg(sigma'), "XYZ");
% translate(target, [30, 10, -20]);
% Target = patch("Faces", target.Faces, "Vertices", target.Vertices);
% chaser = readSurfaceMesh("SmallSat.glb");
% Chaser = patch("Faces", chaser.Faces, "Vertices", chaser.Vertices);
plot3(orbit_m2km(1,:), orbit_m2km(2,:), orbit_m2km(3,:), 'LineWidth', orbit_width);
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


%%
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
targetPos = timeseries(target_state_log(1:3, :)', time_span);
targetPos.Name = pos_name;
targetAtt = timeseries(eul2quat(target_state_log(7:9, :)'), time_span);
targetAtt.Name = att_name;

target_sat = satellite(sc, targetPos, Name='Target_Sat');
pointAt(target_sat, targetAtt);
ScenarioViewer = satelliteScenarioViewer(sc, CameraReferenceFrame='Inertial');
target_sat.Visual3DModel = 'SmallSat.glb';
target_sat.Visual3DModelScale = 0.8;
hide(target_sat.Orbit);
play(sc);
camtarget(ScenarioViewer, target_sat);


