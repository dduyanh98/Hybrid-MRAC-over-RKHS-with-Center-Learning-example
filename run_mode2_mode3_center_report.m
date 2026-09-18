function report = run_mode2_mode3_center_report()
%RUN_MODE2_MODE3_CENTER_REPORT Compare current center modes 2 and 3.
% Mode 2 is the old original-box 2x2/3x3/5x5 center selection mode.

out_dir = fullfile(pwd, 'Results');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
old_skip_plot_build = getenv('PENDULUM_SKIP_PLOT_BUILD');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, ...
    old_keep_workspace, old_skip_plot_build));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '0');
setenv('PENDULUM_SKIP_PLOT_BUILD', '1');

fprintf('Running center mode 2...\n');
mode2 = run_one_mode(2);
fprintf('Finished center mode 2: %d updates.\n', mode2.updates);

fprintf('Running center mode 3...\n');
mode3 = run_one_mode(3);
fprintf('Finished center mode 3: %d updates.\n', mode3.updates);

[switch_table, dwell_table] = center_switching_tables(mode2.center_times, ...
    mode2.center_counts, mode2.Tmax);

writetable(switch_table, fullfile(out_dir, 'mode2_center_switching_properties.csv'));
writetable(dwell_table, fullfile(out_dir, 'mode2_center_dwell_times.csv'));

fprintf('Writing tables and error plot...\n');
setenv('PENDULUM_SKIP_PLOT_BUILD', '0');
fig = plot_error_comparison(mode2, mode3);
savefig(fig, fullfile(out_dir, 'mode2_mode3_error_comparison.fig'));
exportgraphics(fig, fullfile(out_dir, 'mode2_mode3_error_comparison.png'), ...
    'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, 'mode2_mode3_error_comparison.pdf'), ...
    'ContentType', 'vector');

report = struct();
report.mode2 = mode2;
report.mode3 = mode3;
report.switch_table = switch_table;
report.dwell_table = dwell_table;
save(fullfile(out_dir, 'mode2_mode3_center_report.mat'), 'report');

disp('Mode 2 center switching properties:');
disp(switch_table);
disp('Mode 2 center dwell times:');
disp(dwell_table);
fprintf('Saved results to %s\n', out_dir);
end

function result = run_one_mode(kernel_mode)
overrides = struct('architecture_mode', 3, 'kernel_mode', kernel_mode);
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('evalin(''base'', ''Main_sim'')');

T = evalin('base', 'T');
theta = evalin('base', 'theta');
theta_ref = evalin('base', 'theta_ref');
theta_dot = evalin('base', 'theta_dot');
theta_ref_dot = evalin('base', 'theta_ref_dot');
center_times = evalin('base', 'center_times');
Xi_hist = evalin('base', 'Xi_hist');
p = evalin('base', 'p');
updates = evalin('base', 'updates');

[T_unique, idx_last] = unique(T, 'last');
[T_plot, order] = sort(T_unique);
idx_last = idx_last(order);

result = struct();
result.kernel_mode = kernel_mode;
result.T = T_plot;
result.theta_error_deg = rad2deg(theta(idx_last) - theta_ref(idx_last));
result.theta_dot_error_deg_s = rad2deg(theta_dot(idx_last) - theta_ref_dot(idx_last));
result.center_times = center_times(:);
result.center_counts = cellfun(@(Xi_k) size(Xi_k, 1), Xi_hist(:));
result.Tmax = p.Tmax;
result.updates = updates;
end

function [switch_table, dwell_table] = center_switching_tables(center_times, center_counts, Tmax)
transitions = [4 9; 9 4; 9 25; 25 9];
transition_names = {'4_to_9'; '9_to_4'; '9_to_25'; '25_to_9'};

min_dwell_s = nan(size(transitions, 1), 1);
max_dwell_s = nan(size(transitions, 1), 1);
n_switches = zeros(size(transitions, 1), 1);

if numel(center_counts) >= 2
    changed_idx = find(diff(center_counts) ~= 0) + 1;
    for i = 1:size(transitions, 1)
        from_count = transitions(i, 1);
        to_count = transitions(i, 2);
        match_idx = changed_idx(center_counts(changed_idx - 1) == from_count & ...
                                center_counts(changed_idx) == to_count);
        dwell = center_times(match_idx) - center_times(match_idx - 1);
        if ~isempty(dwell)
            min_dwell_s(i) = min(dwell);
            max_dwell_s(i) = max(dwell);
            n_switches(i) = numel(dwell);
        end
    end
end

switch_table = table(transition_names, transitions(:,1), transitions(:,2), ...
    n_switches, min_dwell_s, max_dwell_s, ...
    'VariableNames', {'transition', 'from_centers', 'to_centers', ...
                      'n_switches', 'min_time_before_switch_s', ...
                      'max_time_before_switch_s'});

center_values = [4; 9; 25];
avg_stay_s = nan(numel(center_values), 1);
total_stay_s = zeros(numel(center_values), 1);
n_stays = zeros(numel(center_values), 1);

if ~isempty(center_counts)
    interval_start = center_times(:);
    interval_end = [center_times(2:end); Tmax];
    interval_duration = max(interval_end - interval_start, 0);

    for i = 1:numel(center_values)
        idx = center_counts(:) == center_values(i);
        if any(idx)
            total_stay_s(i) = sum(interval_duration(idx));
            n_stays(i) = nnz(idx);
            avg_stay_s(i) = mean(interval_duration(idx));
        end
    end
end

dwell_table = table(center_values, n_stays, total_stay_s, avg_stay_s, ...
    'VariableNames', {'centers', 'n_update_intervals', ...
                      'total_time_s', 'average_time_stays_s'});
end

function fig = plot_error_comparison(mode2, mode3)
fig = figure('Color', 'w', 'Name', 'Mode 2 vs Mode 3 Error Comparison');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
legend_font_size = 15;

nexttile;
plot(mode2.T, mode2.theta_error_deg, 'b-', 'LineWidth', 1.5); hold on;
plot(mode3.T, mode3.theta_error_deg, 'k--', 'LineWidth', 1.5);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\theta$ [deg]', 'Interpreter', 'latex');
legend({'Center mode 2', 'Center mode 3'}, 'Location', 'best', ...
    'Interpreter', 'latex', 'FontSize', legend_font_size);
title('Error in $\theta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on;

nexttile;
plot(mode2.T, mode2.theta_dot_error_deg_s, 'b-', 'LineWidth', 1.5); hold on;
plot(mode3.T, mode3.theta_dot_error_deg_s, 'k--', 'LineWidth', 1.5);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_{\dot{\theta}}$ [deg/s]', 'Interpreter', 'latex');
legend({'Center mode 2', 'Center mode 3'}, 'Location', 'best', ...
    'Interpreter', 'latex', 'FontSize', legend_font_size);
title('Error in $\dot{\theta}$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
grid on;
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace, old_skip_plot_build)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
setenv('PENDULUM_SKIP_PLOT_BUILD', old_skip_plot_build);
end
