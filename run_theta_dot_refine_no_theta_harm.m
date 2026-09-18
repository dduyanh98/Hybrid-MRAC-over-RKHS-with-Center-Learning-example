function TBL = run_theta_dot_refine_no_theta_harm()
%RUN_THETA_DOT_REFINE_NO_THETA_HARM Reduce theta-dot oscillation safely.
%
% Goal:
%   Keep theta tracking essentially as good as the current theta-priority
%   tuning, then improve theta-dot tracking/oscillation.
%
% Selection:
%   1) feasible candidates satisfy mode-specific theta constraints
%   2) feasible candidates are sorted by theta_dot_l2_5_10
%   3) infeasible candidates are still reported for context

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

mode1_base = [1200 600 300 200 1200 600];
mode3_base = [450 450 340 150 450 450];

mode1_table = [
    mode1_base
    1200 500 300 200 1200 500
    1200 400 300 200 1200 400
    1200 300 300 200 1200 300
    1200 200 300 200 1200 200
    1200 500 260 200 1200 500
    1200 400 260 200 1200 400
    1200 300 260 200 1200 300
    1200 500 220 175 1200 500
    1200 400 220 175 1200 400
    1000 500 300 200 1000 500
    1000 400 260 175 1000 400
    900  500 260 175 900  500
    900  400 260 175 900  400
    700  500 220 150 700  500
];

mode3_table = [
    mode3_base
    450 350 340 150 450 350
    450 300 340 150 450 300
    450 250 340 150 450 250
    450 200 340 150 450 200
    450 350 300 150 450 350
    450 300 300 150 450 300
    450 250 300 150 450 250
    450 350 260 150 450 350
    450 300 260 150 450 300
    500 350 300 150 500 350
    500 300 260 150 500 300
    550 350 300 175 550 350
    600 350 300 175 600 350
    700 500 300 200 700 500
];

specs = {
    struct('kernel_mode', 1, 'gain_table', mode1_table, ...
           'theta_l2_limit', 0.0068, 'max_theta_limit', 0.55)
    struct('kernel_mode', 3, 'gain_table', mode3_table, ...
           'theta_l2_limit', 0.0065, 'max_theta_limit', 0.50)
};

n_tests = size(mode1_table, 1) + size(mode3_table, 1);
results = repmat(empty_result(), n_tests, 1);
best_feasible = struct();

row = 0;
for s = 1:numel(specs)
    spec = specs{s};
    mode_key = sprintf('mode%d', spec.kernel_mode);

    for g = 1:size(spec.gain_table, 1)
        row = row + 1;
        gains = spec.gain_table(g, :);
        fprintf(['\n===== THETA-DOT REFINE %d/%d: mode=%d ' ...
                 'Gx=[%g %g] Gr=%g Gf=%g Gg=[%g %g] =====\n'], ...
            row, n_tests, spec.kernel_mode, gains(1), gains(2), gains(3), gains(4), gains(5), gains(6));

        results(row).kernel_mode = spec.kernel_mode;
        results(row).Gamma_x_theta = gains(1);
        results(row).Gamma_x_theta_dot = gains(2);
        results(row).Gamma_r = gains(3);
        results(row).Gamma_f = gains(4);
        results(row).Gamma_g_theta = gains(5);
        results(row).Gamma_g_theta_dot = gains(6);

        try
            metrics = run_one_candidate(make_overrides(spec.kernel_mode, gains));
            fields = fieldnames(metrics);
            for i = 1:numel(fields)
                field = fields{i};
                results(row).(field) = metrics.(field);
            end

            results(row).theta_feasible = ...
                results(row).theta_l2_5_10 <= spec.theta_l2_limit && ...
                results(row).max_abs_theta_err_5_10 <= spec.max_theta_limit;

            if results(row).theta_feasible && ...
                    (~isfield(best_feasible, mode_key) || is_better_dot(results(row), best_feasible.(mode_key)))
                best_feasible.(mode_key) = results(row);
            end

            fprintf(['thetaL2=%.6f, dotL2=%.6f, dotZC=%g, dotTV=%.4f, ' ...
                     'thetaRMSE=%.4f deg, dotRMSE=%.4f deg/s, maxTheta=%.4f deg, feasible=%d\n'], ...
                results(row).theta_l2_5_10, results(row).theta_dot_l2_5_10, ...
                results(row).theta_dot_zero_crossings_5_10, results(row).theta_dot_total_variation_5_10, ...
                results(row).theta_rmse_5_10, results(row).theta_dot_rmse_5_10, ...
                results(row).max_abs_theta_err_5_10, results(row).theta_feasible);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, ...
    {'kernel_mode', 'theta_feasible', 'theta_dot_l2_5_10', 'theta_l2_5_10'}, ...
    {'ascend', 'descend', 'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'theta_dot_refine_no_theta_harm.mat'), ...
    'results', 'TBL', 'mode1_table', 'mode3_table', 'best_feasible');

disp(TBL(:, {'kernel_mode','theta_feasible','Gamma_x_theta','Gamma_x_theta_dot', ...
    'Gamma_r','Gamma_f','Gamma_g_theta','Gamma_g_theta_dot', ...
    'theta_l2_5_10','theta_dot_l2_5_10','theta_dot_zero_crossings_5_10', ...
    'theta_dot_total_variation_5_10','theta_rmse_5_10','theta_dot_rmse_5_10', ...
    'max_abs_theta_err_5_10','kernel_updates','error'}));
end

function tf = is_better_dot(a, b)
tf = a.theta_dot_l2_5_10 < b.theta_dot_l2_5_10;
if ~tf && abs(a.theta_dot_l2_5_10 - b.theta_dot_l2_5_10) < 1e-12
    tf = a.theta_dot_total_variation_5_10 < b.theta_dot_total_variation_5_10;
end
if ~tf && abs(a.theta_dot_l2_5_10 - b.theta_dot_l2_5_10) < 1e-12 && ...
        abs(a.theta_dot_total_variation_5_10 - b.theta_dot_total_variation_5_10) < 1e-12
    tf = a.theta_l2_5_10 < b.theta_l2_5_10;
end
end

function overrides = make_overrides(kernel_mode, gains)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Gamma_x', diag(gains(1:2)), ...
    'Gamma_r', gains(3), ...
    'Gamma_f', gains(4), ...
    'Gamma_g', diag(gains(5:6)));
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
metrics.theta_dot_zero_crossings_5_10 = sum(dot_err(1:end-1) .* dot_err(2:end) < 0);
metrics.theta_dot_total_variation_5_10 = sum(abs(diff(rad2deg(dot_err))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.kernel_updates = SIM_METRICS.kernel_updates;
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
    'theta_feasible', false, ...
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
    'theta_dot_zero_crossings_5_10', NaN, ...
    'theta_dot_total_variation_5_10', NaN, ...
    'control_rms', NaN, ...
    'control_max_abs', NaN, ...
    'kernel_updates', NaN, ...
    'error', '');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
