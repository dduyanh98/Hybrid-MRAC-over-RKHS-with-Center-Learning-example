function TBL = run_q_tune_theta_dot_tradeoff()
%RUN_Q_TUNE_THETA_DOT_TRADEOFF Tune Q for theta-dot behavior.
%
% Tracks x_ref, not yd.
% Current adaptive gains are left unchanged.
% Q candidates are diagonal. Selection allows modest theta degradation, then
% prioritizes theta-dot L2 and theta-dot total variation.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

kernel_modes = [1 3];
q_table = [
    1.0   1.0
    1.0   2.0
    1.0   5.0
    1.0  10.0
    1.0  15.0
    1.0  20.0
    1.0  30.0
    1.0  50.0
    0.5   5.0
    0.5  10.0
    0.5  20.0
    2.0  10.0
    2.0  20.0
];

n_tests = numel(kernel_modes) * size(q_table, 1);
results = repmat(empty_result(), n_tests, 1);

row = 0;
for m = 1:numel(kernel_modes)
    kernel_mode = kernel_modes(m);
    for q = 1:size(q_table, 1)
        row = row + 1;
        qdiag = q_table(q, :);
        fprintf('\n===== Q TEST %d/%d: mode=%d Q=diag([%g %g]) =====\n', ...
            row, n_tests, kernel_mode, qdiag(1), qdiag(2));

        results(row).kernel_mode = kernel_mode;
        results(row).Q_theta = qdiag(1);
        results(row).Q_theta_dot = qdiag(2);

        try
            metrics = run_one_candidate(make_overrides(kernel_mode, qdiag));
            fields = fieldnames(metrics);
            for i = 1:numel(fields)
                field = fields{i};
                results(row).(field) = metrics.(field);
            end
            fprintf(['thetaL2=%.6f, dotL2=%.6f, dotTV=%.4f, ' ...
                     'thetaRMSE=%.4f deg, dotRMSE=%.4f deg/s, maxTheta=%.4f deg, updates=%g\n'], ...
                results(row).theta_l2_5_10, results(row).theta_dot_l2_5_10, ...
                results(row).theta_dot_total_variation_5_10, results(row).theta_rmse_5_10, ...
                results(row).theta_dot_rmse_5_10, results(row).max_abs_theta_err_5_10, ...
                results(row).kernel_updates);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

results = mark_feasible(results);
TBL = struct2table(results);
TBL = sortrows(TBL, ...
    {'kernel_mode', 'tradeoff_feasible', 'theta_dot_l2_5_10', 'theta_dot_total_variation_5_10'}, ...
    {'ascend', 'descend', 'ascend', 'ascend'});

best_by_mode = best_feasible_by_mode(results, kernel_modes);

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'q_tune_theta_dot_tradeoff.mat'), ...
    'results', 'TBL', 'q_table', 'best_by_mode');

disp(TBL(:, {'kernel_mode','tradeoff_feasible','Q_theta','Q_theta_dot', ...
    'theta_l2_5_10','theta_dot_l2_5_10','theta_dot_total_variation_5_10', ...
    'theta_rmse_5_10','theta_dot_rmse_5_10','max_abs_theta_err_5_10', ...
    'kernel_updates','error'}));
end

function results = mark_feasible(results)
modes = unique([results.kernel_mode]);
for i = 1:numel(modes)
    mode = modes(i);
    idx_mode = find([results.kernel_mode] == mode);
    idx_base = idx_mode([results(idx_mode).Q_theta] == 1 & [results(idx_mode).Q_theta_dot] == 1);
    if isempty(idx_base)
        continue;
    end
    base = results(idx_base(1));
    theta_l2_limit = 1.35 * base.theta_l2_5_10;
    max_theta_limit = max(0.60, 1.35 * base.max_abs_theta_err_5_10);

    for j = idx_mode
        results(j).baseline_theta_l2_5_10 = base.theta_l2_5_10;
        results(j).baseline_theta_dot_l2_5_10 = base.theta_dot_l2_5_10;
        results(j).tradeoff_feasible = ...
            results(j).theta_l2_5_10 <= theta_l2_limit && ...
            results(j).max_abs_theta_err_5_10 <= max_theta_limit;
    end
end
end

function best_by_mode = best_feasible_by_mode(results, kernel_modes)
best_by_mode = struct();
for i = 1:numel(kernel_modes)
    mode = kernel_modes(i);
    idx = find([results.kernel_mode] == mode & [results.tradeoff_feasible]);
    if isempty(idx)
        continue;
    end
    [~, order] = sortrows([[results(idx).theta_dot_l2_5_10]', ...
                           [results(idx).theta_dot_total_variation_5_10]', ...
                           [results(idx).theta_l2_5_10]']);
    best_by_mode.(sprintf('mode%d', mode)) = results(idx(order(1)));
end
end

function overrides = make_overrides(kernel_mode, qdiag)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Q', diag(qdiag));
end

function metrics = run_one_candidate(overrides)
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:,1:2) - X(:,3:4);
idx_late = T >= 5 & T <= 10;
dot_err = e(idx_late, 2);

metrics = struct();
metrics.theta_l2_5_10 = sqrt(trapz(T(idx_late), e(idx_late,1).^2));
metrics.theta_dot_l2_5_10 = sqrt(trapz(T(idx_late), dot_err.^2));
metrics.combined_l2_5_10 = sqrt(trapz(T(idx_late), sum(e(idx_late,:).^2, 2)));
metrics.theta_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,1)).^2));
metrics.theta_dot_rmse_5_10 = sqrt(mean(rad2deg(dot_err).^2));
metrics.max_abs_theta_err_5_10 = max(abs(rad2deg(e(idx_late,1))));
metrics.theta_dot_total_variation_5_10 = sum(abs(diff(rad2deg(dot_err))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.kernel_updates = SIM_METRICS.kernel_updates;
metrics.tradeoff_feasible = false;
metrics.baseline_theta_l2_5_10 = NaN;
metrics.baseline_theta_dot_l2_5_10 = NaN;
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
    'Q_theta', NaN, ...
    'Q_theta_dot', NaN, ...
    'theta_l2_5_10', NaN, ...
    'theta_dot_l2_5_10', NaN, ...
    'combined_l2_5_10', NaN, ...
    'theta_rmse_5_10', NaN, ...
    'theta_dot_rmse_5_10', NaN, ...
    'max_abs_theta_err_5_10', NaN, ...
    'theta_dot_total_variation_5_10', NaN, ...
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'kernel_updates', NaN, ...
    'tradeoff_feasible', false, ...
    'baseline_theta_l2_5_10', NaN, ...
    'baseline_theta_dot_l2_5_10', NaN, ...
    'error', '');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
