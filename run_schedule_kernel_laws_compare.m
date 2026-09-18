function TBL = run_schedule_kernel_laws_compare()
%RUN_SCHEDULE_KERNEL_LAWS_COMPARE Compare mild/current vs aggressive plants.
%
% Candidate set:
%   - plant schedule: aggressive, current
%   - kernel mode: 1, 3
%   - mode 1 center laws: grid sequence 2-3-5 vs 2-4-6
%   - mode 3 center laws: endpoint vs max_tracking_error
%   - relaxed windows: high/low = 0.03/0.008 vs 0.02/0.015
%
% Adaptive gains are not swept here. Each mode uses the gains already defined
% by hybrid_params + apply_kernel_mode_params.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

schedules = {'aggressive', 'current'};
cases = make_cases();
n_tests = numel(schedules) * numel(cases);

results = repmat(empty_result(), n_tests, 1);
best_runs = struct();

row = 0;
for s = 1:numel(schedules)
    schedule = schedules{s};
    fprintf('\n========== SCHEDULE: %s ==========\n', schedule);

    for c = 1:numel(cases)
        row = row + 1;
        test_case = cases(c);
        overrides = make_overrides(schedule, test_case.params);

        results(row).schedule = schedule;
        results(row).case_name = test_case.name;
        results(row).kernel_mode = test_case.params.kernel_mode;
        results(row).kernel_l2_high = test_case.params.kernel_l2_high;
        results(row).kernel_l2_low = test_case.params.kernel_l2_low;

        fprintf('Case %d/%d: %s / %s\n', row, n_tests, schedule, test_case.name);
        try
            run_data = run_one_candidate(overrides);
            metric_fields = fieldnames(run_data.metrics);
            for i = 1:numel(metric_fields)
                field = metric_fields{i};
                results(row).(field) = run_data.metrics.(field);
            end

            if ~isfield(best_runs, schedule) || ...
                    run_data.metrics.l2_5_10 < best_runs.(schedule).metrics.l2_5_10
                best_runs.(schedule) = run_data;
                best_runs.(schedule).schedule = schedule;
                best_runs.(schedule).case_name = test_case.name;
                best_runs.(schedule).overrides = overrides;
            end

            fprintf(['  L2_5_10=%.6f, theta=%.4f deg, thetaDot=%.4f deg/s, ' ...
                     'maxTheta=%.4f deg, updates=%g, centers final/mean/max=%g/%.1f/%g\n'], ...
                results(row).l2_5_10, results(row).theta_rmse_5_10, ...
                results(row).theta_dot_rmse_5_10, results(row).max_abs_theta_err_5_10, ...
                results(row).kernel_updates, results(row).final_centers, ...
                results(row).mean_centers, results(row).max_centers);
        catch ME
            results(row).error = ME.message;
            fprintf('  FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, {'schedule', 'l2_5_10'}, {'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

plot_best_error_comparison(best_runs);
save(fullfile('Results', 'schedule_kernel_laws_compare.mat'), ...
    'results', 'TBL', 'best_runs', 'cases');

disp(TBL(:, {'schedule','case_name','kernel_mode','kernel_l2_high','kernel_l2_low', ...
    'l2_5_10','theta_rmse_5_10','theta_dot_rmse_5_10','max_abs_theta_err_5_10', ...
    'kernel_updates','final_centers','mean_centers','max_centers','error'}));
end

function cases = make_cases()
windows = [
    0.030 0.008
    0.020 0.015
];

case_list = {};
for w = 1:size(windows, 1)
    hi = windows(w, 1);
    lo = windows(w, 2);

    case_list{end+1} = struct('name', sprintf('m1_grid_235_win_%g_%g', hi, lo), ...
        'params', struct('kernel_mode', 1, ...
                         'kernel_mode1_grid_sequence', [2 3 5], ...
                         'kernel_l2_high', hi, ...
                         'kernel_l2_low', lo));

    case_list{end+1} = struct('name', sprintf('m1_grid_246_win_%g_%g', hi, lo), ...
        'params', struct('kernel_mode', 1, ...
                         'kernel_mode1_grid_sequence', [2 4 6], ...
                         'kernel_mode1_max_grid_n', 6, ...
                         'kernel_l2_high', hi, ...
                         'kernel_l2_low', lo));

    case_list{end+1} = struct('name', sprintf('m3_endpoint_win_%g_%g', hi, lo), ...
        'params', struct('kernel_mode', 3, ...
                         'kernel_rep_point', 'endpoint', ...
                         'kernel_l2_high', hi, ...
                         'kernel_l2_low', lo));

    case_list{end+1} = struct('name', sprintf('m3_max_tracking_win_%g_%g', hi, lo), ...
        'params', struct('kernel_mode', 3, ...
                         'kernel_rep_point', 'max_tracking_error', ...
                         'kernel_l2_high', hi, ...
                         'kernel_l2_low', lo));
end

cases = [case_list{:}];
end

function overrides = make_overrides(schedule, params)
overrides = struct( ...
    'plant_schedule', schedule, ...
    'architecture_mode', 3, ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant');

fields = fieldnames(params);
for i = 1:numel(fields)
    field = fields{i};
    overrides.(field) = params.(field);
end
end

function run_data = run_one_candidate(overrides)
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:,1:2) - X(:,3:4);
idx_late = T >= 5 & T <= 10;

metrics = struct();
metrics.l2_5_10 = sqrt(trapz(T(idx_late), sum(e(idx_late,:).^2, 2)));
metrics.l2_full = SIM_METRICS.tracking_l2;
metrics.theta_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,1)).^2));
metrics.theta_dot_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,2)).^2));
metrics.max_abs_theta_err_5_10 = max(abs(rad2deg(e(idx_late,1))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.kernel_updates = SIM_METRICS.kernel_updates;

if exist('Xi_hist', 'var') && ~isempty(Xi_hist)
    center_counts = cellfun(@(Xi) size(Xi, 1), Xi_hist(:));
    metrics.final_centers = center_counts(end);
    metrics.mean_centers = mean(center_counts);
    metrics.max_centers = max(center_counts);
else
    metrics.final_centers = 0;
    metrics.mean_centers = 0;
    metrics.max_centers = 0;
end

run_data = struct();
run_data.T = T;
run_data.e = e;
run_data.metrics = metrics;
end

function plot_best_error_comparison(best_runs)
figure('Color', 'w', 'Visible', 'off');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,1)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,1)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta error [deg]');
title('\theta Tracking Error: Best Aggressive vs Best Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_kernel_laws_theta_error.png'));
savefig(gcf, fullfile('Results', 'schedule_kernel_laws_theta_error.fig'));

figure('Color', 'w', 'Visible', 'off');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,2)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,2)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta dot error [deg/s]');
title('\theta dot Tracking Error: Best Aggressive vs Best Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_kernel_laws_theta_dot_error.png'));
savefig(gcf, fullfile('Results', 'schedule_kernel_laws_theta_dot_error.fig'));
end

function label = best_label(run_data)
label = sprintf('%s best: %s, L2_{5-10}=%.4g', ...
    run_data.schedule, run_data.case_name, run_data.metrics.l2_5_10);
end

function result = empty_result()
result = struct( ...
    'schedule', '', ...
    'case_name', '', ...
    'kernel_mode', NaN, ...
    'kernel_l2_high', NaN, ...
    'kernel_l2_low', NaN, ...
    'l2_5_10', NaN, ...
    'l2_full', NaN, ...
    'theta_rmse_5_10', NaN, ...
    'theta_dot_rmse_5_10', NaN, ...
    'max_abs_theta_err_5_10', NaN, ...
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'kernel_updates', NaN, ...
    'final_centers', NaN, ...
    'mean_centers', NaN, ...
    'max_centers', NaN, ...
    'error', '');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
