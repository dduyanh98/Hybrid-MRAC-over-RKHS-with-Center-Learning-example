function TBL = run_theta_priority_anisotropic_gain_tune()
%RUN_THETA_PRIORITY_ANISOTROPIC_GAIN_TUNE Tune gains with theta-first metric.
%
% Tracks x_ref, not yd. Ranking is:
%   1) theta_l2_5_10
%   2) theta_dot_l2_5_10
%   3) max_abs_theta_err_5_10
%
% Gamma_x and Gamma_g may be anisotropic diagonal matrices:
%   Gamma_x = diag([Gx_theta, Gx_theta_dot])
%   Gamma_g = diag([Gg_theta, Gg_theta_dot])

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

kernel_modes = [1 3];
gain_table = make_gain_table();
n_tests = numel(kernel_modes) * size(gain_table, 1);

results = repmat(empty_result(), n_tests, 1);
best_by_mode = struct();

row = 0;
for m = 1:numel(kernel_modes)
    kernel_mode = kernel_modes(m);
    mode_key = sprintf('mode%d', kernel_mode);

    for g = 1:size(gain_table, 1)
        row = row + 1;
        spec = gain_table(g, :);
        fprintf(['\n===== THETA-FIRST GAIN TEST %d/%d: mode=%d ' ...
                 'Gx=[%g %g] Gr=%g Gf=%g Gg=[%g %g] =====\n'], ...
            row, n_tests, kernel_mode, spec(1), spec(2), spec(3), spec(4), spec(5), spec(6));

        results(row).kernel_mode = kernel_mode;
        results(row).Gamma_x_theta = spec(1);
        results(row).Gamma_x_theta_dot = spec(2);
        results(row).Gamma_r = spec(3);
        results(row).Gamma_f = spec(4);
        results(row).Gamma_g_theta = spec(5);
        results(row).Gamma_g_theta_dot = spec(6);

        try
            metrics = run_one_candidate(make_overrides(kernel_mode, spec));
            fields = fieldnames(metrics);
            for i = 1:numel(fields)
                field = fields{i};
                results(row).(field) = metrics.(field);
            end

            if ~isfield(best_by_mode, mode_key) || is_better(results(row), best_by_mode.(mode_key))
                best_by_mode.(mode_key) = results(row);
            end

            fprintf(['thetaL2=%.6f, thetaDotL2=%.6f, combined=%.6f, ' ...
                     'thetaRMSE=%.4f deg, thetaDotRMSE=%.4f deg/s, maxTheta=%.4f deg, ' ...
                     'updates=%g, centers final/mean/max=%g/%.1f/%g, uRMS=%.4f\n'], ...
                results(row).theta_l2_5_10, results(row).theta_dot_l2_5_10, ...
                results(row).combined_l2_5_10, results(row).theta_rmse_5_10, ...
                results(row).theta_dot_rmse_5_10, results(row).max_abs_theta_err_5_10, ...
                results(row).kernel_updates, results(row).final_centers, ...
                results(row).mean_centers, results(row).max_centers, results(row).control_rms);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, ...
    {'kernel_mode', 'theta_l2_5_10', 'theta_dot_l2_5_10', 'max_abs_theta_err_5_10'}, ...
    {'ascend', 'ascend', 'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'theta_priority_anisotropic_gain_tune.mat'), ...
    'results', 'TBL', 'gain_table', 'best_by_mode');

disp(TBL(:, {'kernel_mode','Gamma_x_theta','Gamma_x_theta_dot','Gamma_r','Gamma_f', ...
    'Gamma_g_theta','Gamma_g_theta_dot','theta_l2_5_10','theta_dot_l2_5_10', ...
    'combined_l2_5_10','theta_rmse_5_10','theta_dot_rmse_5_10', ...
    'max_abs_theta_err_5_10','kernel_updates','final_centers','mean_centers', ...
    'max_centers','control_rms','error'}));
end

function gain_table = make_gain_table()
% [Gx_theta Gx_theta_dot Gamma_r Gamma_f Gg_theta Gg_theta_dot]
gain_table = [
    450  450  340  150  450  450
    500  500  340  175  500  500
    550  550  350  225  550  550
    600  600  450  225  600  600
    700  700  450  275  700  700

    700  250  350  175  700  250
    900  250  350  175  900  250
    1200 250  350  175  1200 250
    900  400  350  200  900  400
    1200 400  350  200  1200 400

    700  200  220  150  700  200
    900  200  220  150  900  200
    1200 200  220  175  1200 200
    900  300  260  175  900  300
    1200 300  260  200  1200 300

    500  200  220  125  900  200
    700  200  220  125  1200 200
    900  200  260  150  1400 200
    1200 250  260  175  1600 250
    1400 300  300  200  1800 300

    700  500  220  150  700  500
    900  500  260  175  900  500
    1200 600 300  200  1200 600
    500  900  260  175  500  900
    700  1200 300  200  700  1200
];
end

function tf = is_better(a, b)
tf = a.theta_l2_5_10 < b.theta_l2_5_10;
if ~tf && abs(a.theta_l2_5_10 - b.theta_l2_5_10) < 1e-12
    tf = a.theta_dot_l2_5_10 < b.theta_dot_l2_5_10;
end
if ~tf && abs(a.theta_l2_5_10 - b.theta_l2_5_10) < 1e-12 && ...
        abs(a.theta_dot_l2_5_10 - b.theta_dot_l2_5_10) < 1e-12
    tf = a.max_abs_theta_err_5_10 < b.max_abs_theta_err_5_10;
end
end

function overrides = make_overrides(kernel_mode, spec)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Gamma_x', diag(spec(1:2)), ...
    'Gamma_r', spec(3), ...
    'Gamma_f', spec(4), ...
    'Gamma_g', diag(spec(5:6)));
end

function metrics = run_one_candidate(overrides)
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:,1:2) - X(:,3:4);
idx_late = T >= 5 & T <= 10;

metrics = struct();
metrics.theta_l2_5_10 = sqrt(trapz(T(idx_late), e(idx_late,1).^2));
metrics.theta_dot_l2_5_10 = sqrt(trapz(T(idx_late), e(idx_late,2).^2));
metrics.combined_l2_5_10 = sqrt(trapz(T(idx_late), sum(e(idx_late,:).^2, 2)));
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
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
    'Gamma_x_theta', NaN, ...
    'Gamma_x_theta_dot', NaN, ...
    'Gamma_r', NaN, ...
    'Gamma_f', NaN, ...
    'Gamma_g_theta', NaN, ...
    'Gamma_g_theta_dot', NaN, ...
    'theta_l2_5_10', NaN, ...
    'theta_dot_l2_5_10', NaN, ...
    'combined_l2_5_10', NaN, ...
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
