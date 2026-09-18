function TBL = run_schedule_tuning_compare()
%RUN_SCHEDULE_TUNING_COMPARE Tune current vs aggressive plant schedules.
%
% Fixed comparison setting:
%   - architecture 3
%   - fixed hybrid budget
%   - constant jump budgets
%   - candidate set varies adaptive gains and kernel center policy
%
% The script saves:
%   Results/schedule_tuning_compare.mat
%   Results/schedule_tuning_theta_error.png
%   Results/schedule_tuning_theta_dot_error.png

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

schedules = {'aggressive', 'current'};
fast_sweep = strcmp(getenv('PENDULUM_FAST_SWEEP'), '1');
kernel_cases = make_kernel_cases(fast_sweep);
gain_table = [
      5    5    5    5
     20   20   10   20
     50   50   40   80
    100  150   50  100
    200  220   75  200
    300  350  150  300
    400  350  150  400
    500  340  175  500
    600  450  225  600
];
if fast_sweep
    gain_table = [
         20   20   10   20
         50   50   40   80
        300  350  150  300
        500  340  175  500
    ];
end

n_tests = numel(schedules) * numel(kernel_cases) * size(gain_table, 1);
results = repmat(empty_result(), n_tests, 1);
best_runs = struct();

row = 0;
for s = 1:numel(schedules)
    schedule = schedules{s};
    fprintf('\n========== TUNING SCHEDULE: %s ==========\n', schedule);

    for k = 1:numel(kernel_cases)
        kernel_case = kernel_cases(k);
        for g = 1:size(gain_table, 1)
            row = row + 1;
            gains = gain_table(g, :);
            overrides = make_overrides(schedule, kernel_case, gains);

            results(row).schedule = schedule;
            results(row).case_name = kernel_case.name;
            results(row).Gamma_x = gains(1);
            results(row).Gamma_r = gains(2);
            results(row).Gamma_f = gains(3);
            results(row).Gamma_g = gains(4);

            fprintf('Test %d/%d: %s / %s / gx=%g gr=%g gf=%g gg=%g\n', ...
                row, n_tests, schedule, kernel_case.name, gains(1), gains(2), gains(3), gains(4));

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
                    best_runs.(schedule).case_name = kernel_case.name;
                    best_runs.(schedule).gains = gains;
                    best_runs.(schedule).overrides = overrides;
                end

                fprintf(['  L2_5_10=%.6f, theta=%.4f deg, thetaDot=%.4f deg/s, ' ...
                         'maxTheta=%.4f deg, updates=%g, centers max=%g, uRMS=%.4f\n'], ...
                    results(row).l2_5_10, results(row).theta_rmse_5_10, ...
                    results(row).theta_dot_rmse_5_10, results(row).max_abs_theta_err_5_10, ...
                    results(row).kernel_updates, results(row).max_centers, results(row).control_rms);
            catch ME
                results(row).error = ME.message;
                fprintf('  FAILED: %s\n', ME.message);
            end
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, {'schedule', 'l2_5_10'}, {'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

plot_best_error_comparison(best_runs);

save(fullfile('Results', 'schedule_tuning_compare.mat'), ...
    'results', 'TBL', 'best_runs', 'kernel_cases', 'gain_table');

disp(TBL(:, {'schedule','case_name','Gamma_x','Gamma_r','Gamma_f','Gamma_g', ...
    'l2_5_10','theta_rmse_5_10','theta_dot_rmse_5_10','max_abs_theta_err_5_10', ...
    'kernel_updates','max_centers','control_rms','error'}));
end

function kernel_cases = make_kernel_cases(fast_sweep)
all_cases = [
    struct('name', 'm1_fast_2_4_6', ...
           'params', struct('kernel_mode', 1, ...
                            'kernel_mode1_grid_sequence', [2 4 6], ...
                            'kernel_mode1_max_grid_n', 6, ...
                            'kernel_max_centers', 40));
    struct('name', 'm1_dense_6_relaxed', ...
           'params', struct('kernel_mode', 1, ...
                            'kernel_mode1_grid_sequence', 6, ...
                            'kernel_mode1_max_grid_n', 6, ...
                            'kernel_l2_high', 0.010, ...
                            'kernel_l2_low', 0.020, ...
                            'kernel_max_centers', 40));
    struct('name', 'm3_endpoint_base', ...
           'params', struct('kernel_mode', 3, ...
                            'kernel_rep_point', 'endpoint'));
    struct('name', 'm3_max_tracking', ...
           'params', struct('kernel_mode', 3, ...
                            'kernel_rep_point', 'max_tracking_error'));
    struct('name', 'm3_max_relaxed_010_020', ...
           'params', struct('kernel_mode', 3, ...
                            'kernel_rep_point', 'max_tracking_error', ...
                            'kernel_l2_high', 0.010, ...
                            'kernel_l2_low', 0.020, ...
                            'max_depth', 7));
    struct('name', 'm3_max_refine_006_012', ...
           'params', struct('kernel_mode', 3, ...
                            'kernel_rep_point', 'max_tracking_error', ...
                            'kernel_l2_high', 0.006, ...
                            'kernel_l2_low', 0.012, ...
                            'max_depth', 8));
];
if fast_sweep
    kernel_cases = all_cases([1 4 5]);
else
    kernel_cases = all_cases;
end
end

function overrides = make_overrides(schedule, kernel_case, gains)
overrides = struct( ...
    'plant_schedule', schedule, ...
    'architecture_mode', 3, ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Gamma_x', gains(1) * eye(2), ...
    'Gamma_r', gains(2), ...
    'Gamma_f', gains(3), ...
    'Gamma_g', gains(4) * eye(2));

fields = fieldnames(kernel_case.params);
for i = 1:numel(fields)
    field = fields{i};
    overrides.(field) = kernel_case.params.(field);
end
end

function run_data = run_one_candidate(overrides)
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:,1:2) - X(:,3:4);
idx_late = T >= 5 & T <= 10;
idx_mid = T >= 2 & T < 5;

metrics = struct();
metrics.l2_5_10 = sqrt(trapz(T(idx_late), sum(e(idx_late,:).^2, 2)));
metrics.l2_2_5 = sqrt(trapz(T(idx_mid), sum(e(idx_mid,:).^2, 2)));
metrics.l2_full = SIM_METRICS.tracking_l2;
metrics.theta_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,1)).^2));
metrics.theta_dot_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,2)).^2));
metrics.max_abs_theta_err_5_10 = max(abs(rad2deg(e(idx_late,1))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.final_V_main = SIM_METRICS.final_V_main;
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
figure('Color', 'w');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,1)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,1)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta error [deg]');
title('\theta Tracking Error: Tuned Aggressive vs Tuned Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_tuning_theta_error.png'));
savefig(gcf, fullfile('Results', 'schedule_tuning_theta_error.fig'));

figure('Color', 'w');
hold on;
plot(best_runs.aggressive.T, rad2deg(best_runs.aggressive.e(:,2)), ...
    'r-', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.aggressive));
plot(best_runs.current.T, rad2deg(best_runs.current.e(:,2)), ...
    'b--', 'LineWidth', 1.6, 'DisplayName', best_label(best_runs.current));
xlabel('Time [s]');
ylabel('\theta dot error [deg/s]');
title('\theta dot Tracking Error: Tuned Aggressive vs Tuned Mild/Current');
legend('show', 'Location', 'best');
grid on;
saveas(gcf, fullfile('Results', 'schedule_tuning_theta_dot_error.png'));
savefig(gcf, fullfile('Results', 'schedule_tuning_theta_dot_error.fig'));
end

function label = best_label(run_data)
label = sprintf('%s best: %s, G=[%g %g %g %g], L2_{5-10}=%.4g', ...
    run_data.schedule, run_data.case_name, run_data.gains(1), ...
    run_data.gains(2), run_data.gains(3), run_data.gains(4), ...
    run_data.metrics.l2_5_10);
end

function result = empty_result()
result = struct( ...
    'schedule', '', ...
    'case_name', '', ...
    'Gamma_x', NaN, ...
    'Gamma_r', NaN, ...
    'Gamma_f', NaN, ...
    'Gamma_g', NaN, ...
    'l2_5_10', NaN, ...
    'l2_2_5', NaN, ...
    'l2_full', NaN, ...
    'theta_rmse_5_10', NaN, ...
    'theta_dot_rmse_5_10', NaN, ...
    'max_abs_theta_err_5_10', NaN, ...
    'final_centers', NaN, ...
    'mean_centers', NaN, ...
    'max_centers', NaN, ...
    'kernel_updates', NaN, ...
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'final_V_main', NaN, ...
    'error', '');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
