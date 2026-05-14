close all; clear; clc;

addpath(genpath(pwd));

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 12); 

load("MC_Results_HOCBF.mat");
hocbf_loggers = mc_loggers;
load("MC_Results_CCBF.mat");
ccbf_loggers = mc_loggers;
simCfg = SimCfg();

N_time = length(simCfg.sim_time);
N_mc = 100;

rho_hocbf_all = zeros(3, N_time, N_mc);
rho_ccbf_all = zeros(3, N_time, N_mc);
V_rho_hocbf_all = zeros(1, N_time, N_mc);
V_vel_hocbf_all = zeros(1, N_time, N_mc);
h_hocbf_all = zeros(1, N_time, N_mc);
V_rho_ccbf_all = zeros(1, N_time, N_mc);
V_vel_ccbf_all = zeros(1, N_time, N_mc);
h_ccbf_all = zeros(1, N_time, N_mc);

for mc_idx = 1:N_mc
    rho_hocbf_all(:,:,mc_idx) = hocbf_loggers{mc_idx}.loggerRelative.state.log(7:9, :);
    rho_ccbf_all(:,:,mc_idx) = ccbf_loggers{mc_idx}.loggerRelative.state.log(7:9, :);
    
    V_rho_hocbf_all(:,:,mc_idx) = hocbf_loggers{mc_idx}.loggerLyapunov.V.log(1,:);
    V_vel_hocbf_all(:,:,mc_idx) = hocbf_loggers{mc_idx}.loggerLyapunov.V.log(2,:);
    h_hocbf_all(:,:,mc_idx) = hocbf_loggers{mc_idx}.loggerBarrier.h.log(1,:);
    
    V_rho_ccbf_all(:,:,mc_idx) = ccbf_loggers{mc_idx}.loggerLyapunov.V.log(1,:);
    V_vel_ccbf_all(:,:,mc_idx) = ccbf_loggers{mc_idx}.loggerLyapunov.V.log(2,:);
    h_ccbf_all(:,:,mc_idx) = ccbf_loggers{mc_idx}.loggerBarrier.h.log(1,:);
end

mean_hocbf = mean(rho_hocbf_all, 3);
std_hocbf = std(rho_hocbf_all, 0, 3);
mean_ccbf = mean(rho_ccbf_all, 3);
std_ccbf = std(rho_ccbf_all, 0, 3);

upper_hocbf = mean_hocbf + 3 * std_hocbf;
lower_hocbf = mean_hocbf - 3 * std_hocbf;
upper_ccbf = mean_ccbf + 3 * std_ccbf;
lower_ccbf = mean_ccbf - 3 * std_ccbf;

mean_V_rho_hocbf = mean(V_rho_hocbf_all, 3);
mean_V_vel_hocbf = mean(V_vel_hocbf_all, 3);

mean_V_rho_ccbf = mean(V_rho_ccbf_all, 3);
mean_V_vel_ccbf = mean(V_vel_ccbf_all, 3);

upper_V_rho_hocbf = max(V_rho_hocbf_all, [], 3);
lower_V_rho_hocbf = min(V_rho_hocbf_all, [], 3);

upper_V_rho_ccbf = max(V_rho_ccbf_all, [], 3);
lower_V_rho_ccbf = min(V_rho_ccbf_all, [], 3);

upper_V_vel_hocbf = max(V_vel_hocbf_all, [], 3);
lower_V_vel_hocbf = min(V_vel_hocbf_all, [], 3);

upper_V_vel_ccbf = max(V_vel_ccbf_all, [], 3);
lower_V_vel_ccbf = min(V_vel_ccbf_all, [], 3);

mean_h_hocbf = mean(h_hocbf_all, 3);
upper_h_hocbf = max(h_hocbf_all, [], 3);
lower_h_hocbf = min(h_hocbf_all, [], 3);

mean_h_ccbf = mean(h_ccbf_all, 3);
upper_h_ccbf = max(h_ccbf_all, [], 3);
lower_h_ccbf = min(h_ccbf_all, [], 3);

line_width = 1.2;
color_ccbf_rgb = [0, 0.4471, 0.7412];
color_hocbf_rgb = [0.8510, 0.3255, 0.0980];
alpha_val = 0.4;
fig_size = [600, 675];
t = simCfg.sim_time;
t_fill = [t, fliplr(t)];

rhoFig = figure('Theme', 'light');
rhoFig.Position(3:4) = fig_size;
t_rho = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_rho = {'$\rho_x$ (m)', '$\rho_y$ (m)', '$\rho_z$ (m)'};
rho_targets = [1, 0, 0];

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    
    y_fill_hocbf = [lower_hocbf(ax_idx, :), fliplr(upper_hocbf(ax_idx, :))];
    fill(ax, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    
    y_fill_ccbf = [lower_ccbf(ax_idx, :), fliplr(upper_ccbf(ax_idx, :))];
    fill(ax, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    
    p_hocbf_sample = plot(ax, t, mean_hocbf(ax_idx, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
    p_ccbf_sample = plot(ax, t, mean_ccbf(ax_idx, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
    
    p_target = yline(ax, rho_targets(ax_idx), 'k--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_rho{ax_idx}, 'Interpreter', 'latex');
    
    if ax_idx == 1
        legend(ax, [p_ccbf_sample, p_hocbf_sample, p_target], ...
               {'CCBF', 'HOCBF', 'Target'}, ...
               'Location', 'northeast', 'NumColumns', 3);
    end
    if ax_idx == 3
        xlabel(ax, 'Time (s)'); 
    end
end

rhoZoomFig = figure('Theme', 'light');
rhoZoomFig.Position(3:4) = fig_size;
t_zoom_rho = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylims = {[0, 5], [-2.5, 2.5], [-2.5, 2.5]};

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    
    y_fill_hocbf = [lower_hocbf(ax_idx, :), fliplr(upper_hocbf(ax_idx, :))];
    fill(ax, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    
    y_fill_ccbf = [lower_ccbf(ax_idx, :), fliplr(upper_ccbf(ax_idx, :))];
    fill(ax, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
    
    p_hocbf_sample = plot(ax, t, mean_hocbf(ax_idx, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
    p_ccbf_sample = plot(ax, t, mean_ccbf(ax_idx, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
    
    p_target = yline(ax, rho_targets(ax_idx), 'k--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_rho{ax_idx}, 'Interpreter', 'latex');
    ylim(ylims{ax_idx});
    
    if ax_idx == 1
        legend(ax, [p_ccbf_sample, p_hocbf_sample, p_target], ...
               {'CCBF', 'HOCBF', 'Target'}, ...
               'Location', 'northeast', 'NumColumns', 3);
    end
    if ax_idx == 3
        xlabel(ax, 'Time (s)'); 
    end
end

lyapFig = figure('Theme', 'light');
lyapFig.Position(3:4) = fig_size;
t_lyap = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
y_fill_hocbf = [lower_V_rho_hocbf(1, :), fliplr(upper_V_rho_hocbf(1, :))];
fill(ax, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
y_fill_ccbf = [lower_V_rho_ccbf(1, :), fliplr(upper_V_rho_ccbf(1, :))];
fill(ax, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
p_hocbf = plot(ax, t, mean_V_rho_hocbf(1, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
p_ccbf = plot(ax, t, mean_V_rho_ccbf(1, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
yline(ax, 0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.8);
ylabel(ax, '$V_{\rho}$', 'Interpreter', 'latex');
legend(ax, [p_ccbf, p_hocbf], {'CCBF', 'HOCBF'}, 'Location', 'northeast', 'NumColumns', 2);

ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
y_fill_hocbf = [lower_V_vel_hocbf(1, :), fliplr(upper_V_vel_hocbf(1, :))];
fill(ax, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
y_fill_ccbf = [lower_V_vel_ccbf(1, :), fliplr(upper_V_vel_ccbf(1, :))];
fill(ax, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
plot(ax, t, mean_V_vel_hocbf(1, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
plot(ax, t, mean_V_vel_ccbf(1, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
yline(ax, 0, 'k--', 'LineWidth', 1.5, 'Alpha', 0.8);
ylabel(ax, '$V_v$', 'Interpreter', 'latex');
xlabel(ax, 'Time (s)');

barrierFig = figure('Theme', 'light');
barrierFig.Position(3:4) = fig_size;
t_barrier = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_h1 = nexttile([2, 1]); hold(ax_h1, 'on'); grid(ax_h1, 'on');
y_fill_hocbf = [lower_h_hocbf(1, :), fliplr(upper_h_hocbf(1, :))];
fill(ax_h1, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
y_fill_ccbf = [lower_h_ccbf(1, :), fliplr(upper_h_ccbf(1, :))];
fill(ax_h1, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
p_hocbf = plot(ax_h1, t, mean_h_hocbf(1, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
p_ccbf = plot(ax_h1, t, mean_h_ccbf(1, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
yline(ax_h1, 0, 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
ylabel(ax_h1, '$h(\mathbf{x})$', 'Interpreter', 'latex');
legend(ax_h1, [p_ccbf, p_hocbf], {'CCBF', 'HOCBF'}, 'Location', 'northeast', 'NumColumns', 2);

ax_h2 = nexttile; hold(ax_h2, 'on'); grid(ax_h2, 'on');
fill(ax_h2, t_fill, y_fill_hocbf, color_hocbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
fill(ax_h2, t_fill, y_fill_ccbf, color_ccbf_rgb, 'FaceAlpha', alpha_val, 'EdgeColor', 'none');
plot(ax_h2, t, mean_h_hocbf(1, :), 'Color', color_hocbf_rgb, 'LineWidth', line_width);
plot(ax_h2, t, mean_h_ccbf(1, :), 'Color', color_ccbf_rgb, 'LineWidth', line_width);
yline(ax_h2, 0, 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
ylim(ax_h2, [-1, 5]); 
ylabel(ax_h2, '$h(\mathbf{x})$ (Zoom)', 'Interpreter', 'latex');
xlabel(ax_h2, 'Time (s)');


dt = simCfg.sim_time(2) - simCfg.sim_time(1);
N_mc = length(hocbf_loggers);

path_length_hocbf = zeros(N_mc, 1);
path_length_ccbf = zeros(N_mc, 1);
effort_hocbf = zeros(N_mc, 1);
effort_ccbf = zeros(N_mc, 1);

for mc_idx = 1:N_mc
    rho_h = hocbf_loggers{mc_idx}.loggerRelative.state.log(7:9, :);
    F_h = hocbf_loggers{mc_idx}.loggerControl.force.log;
    
    rho_c = ccbf_loggers{mc_idx}.loggerRelative.state.log(7:9, :);
    F_c = ccbf_loggers{mc_idx}.loggerControl.force.log;
    
    path_length_hocbf(mc_idx) = sum(vecnorm(diff(rho_h, 1, 2), 2, 1));
    effort_hocbf(mc_idx) = sum(vecnorm(F_h, 2, 1).^2) * dt;
    
    path_length_ccbf(mc_idx) = sum(vecnorm(diff(rho_c, 1, 2), 2, 1));
    effort_ccbf(mc_idx) = sum(vecnorm(F_c, 2, 1).^2) * dt;
end

fprintf('====================================================\n');
fprintf('[Monte-Carlo Simulation Results (%d runs)]\n', N_mc);
fprintf('====================================================\n');
fprintf('Path Length (m):\n');
fprintf('  HOCBF : %.2f +/- %.2f\n', mean(path_length_hocbf), std(path_length_hocbf));
fprintf('  CCBF  : %.2f +/- %.2f\n', mean(path_length_ccbf), std(path_length_ccbf));
fprintf('----------------------------------------------------\n');
fprintf('Control Effort (N^2 s):\n');
fprintf('  HOCBF : %.2f +/- %.2f\n', mean(effort_hocbf), std(effort_hocbf));
fprintf('  CCBF  : %.2f +/- %.2f\n', mean(effort_ccbf), std(effort_ccbf));
fprintf('====================================================\n\n');

color_ccbf_rgb = [0, 0.4471, 0.7412];
color_hocbf_rgb = [0.8510, 0.3255, 0.0980];


boxFig = figure('Theme', 'light');
boxFig.Position(3:4) = [1000, 500];
t_box = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
b1_c = boxchart(ax1, ones(N_mc,1), path_length_ccbf, 'BoxFaceColor', color_ccbf_rgb, 'MarkerColor', color_ccbf_rgb);
b1_h = boxchart(ax1, 2*ones(N_mc,1), path_length_hocbf, 'BoxFaceColor', color_hocbf_rgb, 'MarkerColor', color_hocbf_rgb);
xticks(ax1, [1, 2]);
xticklabels(ax1, {'CCBF', 'HOCBF'});
ylabel(ax1, 'Total Path Length (m)');
xlim(ax1, [0.5, 2.5]);

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
b2_c = boxchart(ax2, ones(N_mc,1), effort_ccbf, 'BoxFaceColor', color_ccbf_rgb, 'MarkerColor', color_ccbf_rgb);
b2_h = boxchart(ax2, 2*ones(N_mc,1), effort_hocbf, 'BoxFaceColor', color_hocbf_rgb, 'MarkerColor', color_hocbf_rgb);
xticks(ax2, [1, 2]);
xticklabels(ax2, {'CCBF', 'HOCBF'});
ylabel(ax2, 'Control Effort $\int \|F\|_2^2 dt$ ($N^2 \cdot s$)');
xlim(ax2, [0.5, 2.5]);


save_dir = 'assets';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

figs_to_save = { ...
    rhoFig,     '01_MC_Relative_Position'; ...
    rhoZoomFig, '02_MC_Relative_Position_Zoom'; ...
    lyapFig,    '03_MC_Lyapunov_Function'; ...
    barrierFig, '04_MC_Barrier_Function'; ...
    boxFig,     '05_MC_Box_Plot_Path_and_Force'
};

for idx = 1:size(figs_to_save, 1)
    fig_handle = figs_to_save{idx, 1};
    file_name  = figs_to_save{idx, 2};
    png_path = fullfile(save_dir, [file_name, '.png']);
    exportgraphics(fig_handle, png_path, 'Resolution', 300);
    fprintf('Saved: %s.png\n', file_name);
end
%%
print_box_stats('Path Length (CCBF)', path_length_ccbf);
print_box_stats('Path Length (HOCBF)', path_length_hocbf);
print_box_stats('Control Effort (CCBF)', effort_ccbf);
print_box_stats('Control Effort (HOCBF)', effort_hocbf);


function print_box_stats(data_name, data)
    q1 = prctile(data, 25);
    q2 = median(data);
    q3 = prctile(data, 75);
    iqr_val = q3 - q1;

    % 수염(Whisker)의 경계값 (Q1 - 1.5*IQR 및 Q3 + 1.5*IQR 내부의 실제 데이터)
    lower_whisker = min(data(data >= q1 - 1.5 * iqr_val));
    upper_whisker = max(data(data <= q3 + 1.5 * iqr_val));

    % 이상치(Outliers)
    outliers = data(data < q1 - 1.5 * iqr_val | data > q3 + 1.5 * iqr_val);

    fprintf('[%s]\n', data_name);
    fprintf('  - Upper Whisker (최댓값 경계): %.2f\n', upper_whisker);
    fprintf('  - Q3 (상위 25%%)            : %.2f\n', q3);
    fprintf('  - Median (중앙값)          : %.2f\n', q2);
    fprintf('  - Q1 (하위 25%%)            : %.2f\n', q1);
    fprintf('  - Lower Whisker (최솟값 경계): %.2f\n', lower_whisker);
    fprintf('  - Outliers (이상치 개수)   : %d개\n', length(outliers));
    fprintf('----------------------------------------------------\n');
end