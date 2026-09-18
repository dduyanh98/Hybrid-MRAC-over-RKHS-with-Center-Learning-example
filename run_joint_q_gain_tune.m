function TBL = run_joint_q_gain_tune()
%RUN_JOINT_Q_GAIN_TUNE Tune Q and adaptive gains together.
%
% Tracks x_ref, not yd. Selection prioritizes theta-dot improvement while
% keeping theta in a reasonable band.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

mode_specs = {
    struct('kernel_mode', 1, 'candidates', mode1_candidates(), ...
           'theta_l2_limit', 0.020, 'max_theta_limit', 0.90)
    struct('kernel_mode', 3, 'candidates', mode3_candidates(), ...
           'theta_l2_limit', 0.020, 'max_theta_limit', 0.90)
};

n_tests = 0;
for i = 1:numel(mode_specs)
    n_tests = n_tests + size(mode_specs{i}.candidates, 1);
end

results = repmat(empty_result(), n_tests, 1);
best_feasible = struct();

row = 0;
for s = 1:numel(mode_specs)
    spec = mode_specs{s};
    mode_key = sprintf('mode%d', spec.kernel_mode);

    for c = 1:size(spec.candidates, 1)
        row = row + 1;
        cand = spec.candidates(c, :);
        qdiag = cand(1:2);
        gains = cand(3:8);

        fprintf(['\n===== JOINT Q/GAIN TEST %d/%d: mode=%d Q=[%g %g] ' ...
                 'Gx=[%g %g] Gr=%g Gf=%g Gg=[%g %g] =====\n'], ...
            row, n_tests, spec.kernel_mode, qdiag(1), qdiag(2), ...
            gains(1), gains(2), gains(3), gains(4), gains(5), gains(6));

        results(row).kernel_mode = spec.kernel_mode;
        results(row).Q_theta = qdiag(1);
        results(row).Q_theta_dot = qdiag(2);
        results(row).Gamma_x_theta = gains(1);
        results(row).Gamma_x_theta_dot = gains(2);
        results(row).Gamma_r = gains(3);
        results(row).Gamma_f = gains(4);
        results(row).Gamma_g_theta = gains(5);
        results(row).Gamma_g_theta_dot = gains(6);

        try
            metrics = run_one_candidate(make_overrides(spec.kernel_mode, qdiag, gains));
            fields = fieldnames(metrics);
            for i = 1:numel(fields)
                field = fields{i};
                results(row).(field) = metrics.(field);
            end

            results(row).tradeoff_feasible = ...
                results(row).theta_l2_5_10 <= spec.theta_l2_limit && ...
                results(row).max_abs_theta_err_5_10 <= spec.max_theta_limit;

            if results(row).tradeoff_feasible && ...
                    (~isfield(best_feasible, mode_key) || is_better(results(row), best_feasible.(mode_key)))
                best_feasible.(mode_key) = results(row);
            end

            fprintf(['thetaL2=%.6f, dotL2=%.6f, dotTV=%.2f, ' ...
                     'thetaRMSE=%.4f deg, dotRMSE=%.4f deg/s, maxTheta=%.4f deg, feasible=%d\n'], ...
                results(row).theta_l2_5_10, results(row).theta_dot_l2_5_10, ...
                results(row).theta_dot_total_variation_5_10, results(row).theta_rmse_5_10, ...
                results(row).theta_dot_rmse_5_10, results(row).max_abs_theta_err_5_10, ...
                results(row).tradeoff_feasible);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, ...
    {'kernel_mode', 'tradeoff_feasible', 'theta_dot_l2_5_10', 'theta_l2_5_10'}, ...
    {'ascend', 'descend', 'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'joint_q_gain_tune.mat'), ...
    'results', 'TBL', 'mode_specs', 'best_feasible');

disp(TBL(:, {'kernel_mode','tradeoff_feasible','Q_theta','Q_theta_dot', ...
    'Gamma_x_theta','Gamma_x_theta_dot','Gamma_r','Gamma_f', ...
    'Gamma_g_theta','Gamma_g_theta_dot','theta_l2_5_10','theta_dot_l2_5_10', ...
    'theta_dot_total_variation_5_10','theta_rmse_5_10','theta_dot_rmse_5_10', ...
    'max_abs_theta_err_5_10','kernel_updates','error'}));
end

function C = mode1_candidates()
% [Qtheta Qdot GxTheta GxDot Gr Gf GgTheta GgDot]
C = [
    1.0 1.0 1000 500 300 200 1000 500
    1.0 2.0 1000 500 300 200 1000 500
    1.0 2.0 1200 500 300 200 1200 500
    1.0 2.0 1400 500 260 200 1400 500
    1.0 2.0 1600 600 260 225 1600 600
    1.0 2.0 1800 600 220 225 1800 600
    1.0 3.0 1400 500 260 200 1400 500
    1.0 3.0 1800 600 220 225 1800 600
    1.0 5.0 1800 600 220 225 1800 600
    1.0 5.0 2200 700 220 250 2200 700
    2.0 5.0 1400 500 260 200 1400 500
    2.0 10.0 1600 600 220 225 1600 600
    2.0 20.0 1600 600 220 225 1600 600
    0.5 5.0 1800 600 220 225 1800 600
    0.5 10.0 2200 700 220 250 2200 700
    1.0 10.0 2200 700 180 250 2200 700
    1.0 20.0 2200 700 180 250 2200 700
    2.0 10.0 2200 700 180 250 2200 700
];
end

function C = mode3_candidates()
% [Qtheta Qdot GxTheta GxDot Gr Gf GgTheta GgDot]
C = [
    1.0 1.0 450 300 340 150 450 300
    1.0 2.0 450 300 340 150 450 300
    1.0 3.0 450 300 340 150 450 300
    1.0 5.0 450 300 340 150 450 300
    1.0 5.0 500 300 300 150 500 300
    1.0 5.0 600 350 300 175 600 350
    1.0 5.0 700 400 260 175 700 400
    1.0 5.0 900 500 260 200 900 500
    1.0 10.0 700 400 260 175 700 400
    1.0 10.0 900 500 220 200 900 500
    2.0 10.0 700 400 260 175 700 400
    2.0 20.0 700 400 220 175 700 400
    0.5 5.0 900 500 260 200 900 500
    0.5 10.0 900 500 220 200 900 500
    1.0 15.0 900 500 220 200 900 500
    1.0 20.0 900 500 220 200 900 500
    2.0 20.0 900 500 220 200 900 500
    1.0 50.0 900 500 220 200 900 500
];
end

function tf = is_better(a, b)
tf = a.theta_dot_l2_5_10 < b.theta_dot_l2_5_10;
if ~tf && abs(a.theta_dot_l2_5_10 - b.theta_dot_l2_5_10) < 1e-12
    tf = a.theta_dot_total_variation_5_10 < b.theta_dot_total_variation_5_10;
end
if ~tf && abs(a.theta_dot_l2_5_10 - b.theta_dot_l2_5_10) < 1e-12 && ...
        abs(a.theta_dot_total_variation_5_10 - b.theta_dot_total_variation_5_10) < 1e-12
    tf = a.theta_l2_5_10 < b.theta_l2_5_10;
end
end

function overrides = make_overrides(kernel_mode, qdiag, gains)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'Q', diag(qdiag), ...
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
metrics.theta_dot_total_variation_5_10 = sum(abs(diff(rad2deg(dot_err))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.kernel_updates = SIM_METRICS.kernel_updates;
metrics.tradeoff_feasible = false;
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
    'tradeoff_feasible', false, ...
    'Q_theta', NaN, ...
    'Q_theta_dot', NaN, ...
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
