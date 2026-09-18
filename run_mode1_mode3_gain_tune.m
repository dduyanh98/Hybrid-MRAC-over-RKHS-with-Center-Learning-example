function TBL = run_mode1_mode3_gain_tune()
%RUN_MODE1_MODE3_GAIN_TUNE Tune adaptive gains for kernel modes 1 and 3.
%
% Fixed setup follows the current project settings:
%   - architecture_mode = 3
%   - plant_schedule = current
%   - hybrid_budget_mode = fixed
%   - fixed_budget_law = constant
%   - kernel mode parameters from hybrid_params/apply_kernel_mode_params
%
% Ranking metric is late tracking L2 over 5-10 seconds.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

kernel_modes = [1 3];
gain_table = [
      5    5    4    8
     20   20   10   20
     50   50   40   80
     50  100   40   80
     80  120   40   80
    100  150   50  100
    150  150   50  150
    200  150   75  200
    200  220   75  200
    250  260  100  250
    300  320  125  300
    300  350  150  300
    350  340  125  350
    400  320  150  400
    400  350  150  400
    450  340  150  450
    500  320  200  500
    500  340  175  500
    550  350  225  550
    600  450  225  600
];

n_tests = numel(kernel_modes) * size(gain_table, 1);
results = repmat(empty_result(), n_tests, 1);
best_by_mode = struct();

row = 0;
for m = 1:numel(kernel_modes)
    kernel_mode = kernel_modes(m);
    mode_key = sprintf('mode%d', kernel_mode);

    for g = 1:size(gain_table, 1)
        row = row + 1;
        gains = gain_table(g, :);
        fprintf('\n===== GAIN TEST %d/%d: kernel_mode=%d gx=%g gr=%g gf=%g gg=%g =====\n', ...
            row, n_tests, kernel_mode, gains(1), gains(2), gains(3), gains(4));

        results(row).kernel_mode = kernel_mode;
        results(row).Gamma_x = gains(1);
        results(row).Gamma_r = gains(2);
        results(row).Gamma_f = gains(3);
        results(row).Gamma_g = gains(4);

        try
            run_data = run_one_candidate(make_overrides(kernel_mode, gains));
            metric_fields = fieldnames(run_data.metrics);
            for i = 1:numel(metric_fields)
                field = metric_fields{i};
                results(row).(field) = run_data.metrics.(field);
            end

            if ~isfield(best_by_mode, mode_key) || ...
                    run_data.metrics.l2_5_10 < best_by_mode.(mode_key).metrics.l2_5_10
                best_by_mode.(mode_key) = run_data;
                best_by_mode.(mode_key).kernel_mode = kernel_mode;
                best_by_mode.(mode_key).gains = gains;
            end

            fprintf(['L2_5_10=%.6f, L2_2_5=%.6f, full=%.6f, ' ...
                     'theta=%.4f deg, thetaDot=%.4f deg/s, maxTheta=%.4f deg, ' ...
                     'updates=%g, centers final/mean/max=%g/%.1f/%g, uRMS=%.4f, uMax=%.4f\n'], ...
                results(row).l2_5_10, results(row).l2_2_5, results(row).l2_full, ...
                results(row).theta_rmse_5_10, results(row).theta_dot_rmse_5_10, ...
                results(row).max_abs_theta_err_5_10, results(row).kernel_updates, ...
                results(row).final_centers, results(row).mean_centers, results(row).max_centers, ...
                results(row).control_rms, results(row).control_max_abs);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, {'kernel_mode', 'l2_5_10'}, {'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'mode1_mode3_gain_tune.mat'), ...
    'results', 'TBL', 'gain_table', 'best_by_mode');

disp(TBL(:, {'kernel_mode','Gamma_x','Gamma_r','Gamma_f','Gamma_g', ...
    'l2_5_10','l2_2_5','l2_full','theta_rmse_5_10','theta_dot_rmse_5_10', ...
    'max_abs_theta_err_5_10','kernel_updates','final_centers','mean_centers', ...
    'max_centers','control_rms','control_max_abs','error'}));
end

function overrides = make_overrides(kernel_mode, gains)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Gamma_x', gains(1) * eye(2), ...
    'Gamma_r', gains(2), ...
    'Gamma_f', gains(3), ...
    'Gamma_g', gains(4) * eye(2));
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

run_data = struct('T', T, 'e', e, 'metrics', metrics);
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
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
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'final_V_main', NaN, ...
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
