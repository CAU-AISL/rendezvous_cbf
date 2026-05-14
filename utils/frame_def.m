clear; close all; clc;

[x, y, z] = sphere(50);
Re = 6.371;
x = x*Re;
y = y*Re;
z = z*Re;
load('topo.mat', 'topo'); 

meshData = readSurfaceMesh('SmallSat.glb');
V = double(meshData.Vertices);
F = double(meshData.Faces);

V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';

scale_factor = 0.8; 
V_base = V * scale_factor;

target_rot = eul2rotm(deg2rad([30, 0, 10]));
target_pos = [10, 0, 3];
targetV = V_base * target_rot + target_pos;
targetF = F;

chaser_rot = eul2rotm(deg2rad([50, 10, -10]));
chaser_pos = [15, 0, 0];
chaserV = V_base * chaser_rot + chaser_pos;
chaserF = F;

FrameFig = figure();
FrameFig.Theme = 'light';
FrameFig.Position = [300, 100, 1631, 1178];
hold on; grid on;

s = surface(x, y, z);
s.FaceColor = 'texturemap';    
s.CData = topo;                
s.EdgeColor = 'none';          
s.FaceAlpha = 0.3;

axis equal;
view([25, 22]);
colormap([zeros(64,1), zeros(64,1), linspace(0.5,1,64)';...
              linspace(0.5,0.8,64)', linspace(1,0.7,64)', zeros(64,1)]);

trisurf(targetF, targetV(:,1), targetV(:,2), targetV(:,3), ...
      'FaceColor', [0.8, 0.8, 0], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');

trisurf(chaserF, chaserV(:,1), chaserV(:,2), chaserV(:,3), ...
      'FaceColor', [0, 0.8, 0.8], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');

camlight('headlight');

I3 = eye(3);
L_E = 8;
labels_E = {'$\hat{i}_i$', '$\hat{j}_i$', '$\hat{k}_i$'};
for i = 1:3
    dir_E = I3(i, :) * L_E;
    dir_E = (eul2rotm(deg2rad([-90, 0, 0]))*dir_E')';
    quiver3(0, 0, 0, dir_E(1), dir_E(2), dir_E(3), ...
        'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 0.3);
end

L_T = 2.5; 
labels_T = {'$\hat{i}_t$', '$\hat{j}_t$', '$\hat{k}_t$'};
for i = 1:3
    dir_T = I3(i, :) * target_rot * L_T; 
    
    quiver3(target_pos(1), target_pos(2), target_pos(3), ...
        dir_T(1), dir_T(2), dir_T(3), ...
        'Color', 'b', 'LineWidth', 2, 'MaxHeadSize', 0.4);
end

L_C = 2.5;
labels_C = {'$\hat{i}_c$', '$\hat{j}_c$', '$\hat{k}_c$'};
for i = 1:3
    dir_C = I3(i, :) * chaser_rot * L_C;
    
    quiver3(chaser_pos(1), chaser_pos(2), chaser_pos(3), ...
        dir_C(1), dir_C(2), dir_C(3), ...
        'Color', 'r', 'LineWidth', 2, 'MaxHeadSize', 0.4);
end
rho_vec = chaser_pos - target_pos;
quiver3(target_pos(1), target_pos(2), target_pos(3), ...
    rho_vec(1), rho_vec(2), rho_vec(3), ...
    'Color', '#00AA00', 'LineWidth', 2.5, 'MaxHeadSize', 0.15, 'LineStyle', '-');

axis equal;
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);

assetsDir = fullfile(pwd, 'assets');
if ~exist(assetsDir, 'dir')
    mkdir(assetsDir);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = fullfile(assetsDir, ['figure_' timestamp '.png']);

set(FrameFig, 'Renderer', 'painters'); 
origUnits = FrameFig.Units;
FrameFig.Units = 'pixels';
figPos = FrameFig.Position;
FrameFig.PaperUnits = 'points';
FrameFig.PaperPosition = [0 0 figPos(3) figPos(4)];
FrameFig.PaperSize = [figPos(3) figPos(4)];
FrameFig.Units = origUnits;

print(FrameFig, filename, '-dpng', '-r300');
fprintf('Saved high-resolution figure to: %s\n', filename);