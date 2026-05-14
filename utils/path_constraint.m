close all; clear; clc;

save_dir = 'assets/';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

offset = 1;
x_min = 1;
x_max = 10 + offset;
x = linspace(x_min, x_max, 100);
theta = linspace(0, 2*pi, 50);
a_h = 0.1;
delta_h = 1;
[X, Theta] = meshgrid(x, theta);
R = sqrt(a_h*(X - delta_h).^3);
Y = R .* cos(Theta);
Z = R .* sin(Theta);

rtPlot = figure();
rtPlot.Theme = 'light';
rtPlot.Position(3:4) = [980, 837];

% [추가] 학술 논문용 폰트/수식 렌더링 설정
set(rtPlot.CurrentAxes, 'TickLabelInterpreter', 'latex');

hold on; grid on;

% Cone 그리기 (투명도 포함)
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'r');

% 위성 메쉬 불러오기
meshData = readSurfaceMesh('SmallSat.glb');
V = double(meshData.Vertices);
F = double(meshData.Faces);
V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';
scale_factor = 0.8; 
V = V * scale_factor;

% 위성 그리기
trisurf(F, V(:,1), V(:,2), V(:,3), ...
      'FaceColor', [0.9, 0.9, 0.9], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');
camlight('headlight');
view([32, 2]);

% 축 라벨 (LaTeX 인터프리터 적용)
xlabel('$x_t$ (m)', 'Interpreter', 'latex'); 
ylabel('$y_t$ (m)', 'Interpreter', 'latex'); 
zlabel('$z_t$ (m)', 'Interpreter', 'latex');

axis equal;
xlim([-2, 10]); ylim([-10, 10]); 

% --- 그림 내보내기 ---
pdf_path = fullfile(save_dir, 'rt_plot.pdf');
png_path = fullfile(save_dir, 'rt_plot.png');

% PDF 저장 (투명도/3D로 인해 자동으로 최적화된 혼합 방식으로 저장됨)
exportgraphics(rtPlot, pdf_path, 'ContentType', 'vector');
fprintf('Saved PDF: %s\n', pdf_path);

% PNG 저장 (300 DPI 고해상도 - PDF가 무거울 경우 대비용)
exportgraphics(rtPlot, png_path, 'Resolution', 300);
fprintf('Saved PNG: %s\n', png_path);