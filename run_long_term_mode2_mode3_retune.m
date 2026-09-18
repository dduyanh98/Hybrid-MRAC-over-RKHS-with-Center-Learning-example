function TBL = run_long_term_mode2_mode3_retune()
%RUN_LONG_TERM_MODE2_MODE3_RETUNE Retune aggressive-plant post-5s L2.
%
% Candidate columns:
% [Atheta Adot Qtheta Qdot GxTheta GxDot Gr Gf GgTheta GgDot hi lo]

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

mode_specs = {
    struct('kernel_mode', 2, 'candidates', mode2_candidates())
    struct('kernel_mode', 3, 'candidates', mode3_candidates())
};

n_tests = 0;
for i = 1:numel(mode_specs)
    n_tests = n_tests + size(mode_specs{i}.candidates, 1);
end

results = repmat(empty_result(), n_tests, 1);
best_by_mode = struct();
row = 0;

for s = 1:numel(mode_specs)
    spec = mode_specs{s};
    mode_key = sprintf('mode%d', spec.kernel_mode);

    for c = 1:size(spec.candidates, 1)
        row = row + 1;
        cand = spec.candidates(c, :);

        fprintf(['\n===== LONG-TAIL RETUNE %d/%d: mode=%d ' ...
                 'A=[-%g -%g] Q=[%g %g] Gx=[%g %g] Gr=%g Gf=%g Gg=[%g %g] hi/lo=[%g %g] =====\n'], ...
            row, n_tests, spec.kernel_mode, cand(1), cand(2), cand(3), cand(4), ...
            cand(5), cand(6), cand(7), cand(8), cand(9), cand(10), cand(11), cand(12));

        results(row).kernel_mode = spec.kernel_mode;
        results(row).A_theta = cand(1);
        results(row).A_theta_dot = cand(2);
        results(row).Q_theta = cand(3);
        results(row).Q_theta_dot = cand(4);
        results(row).Gamma_x_theta = cand(5);
        results(row).Gamma_x_theta_dot = cand(6);
        results(row).Gamma_r = cand(7);
        results(row).Gamma_f = cand(8);
        results(row).Gamma_g_theta = cand(9);
        results(row).Gamma_g_theta_dot = cand(10);
        results(row).kernel_l2_high = cand(11);
        results(row).kernel_l2_low = cand(12);

        try
            metrics = run_one_candidate(make_overrides(spec.kernel_mode, cand));
            fields = fieldnames(metrics);
            for i = 1:numel(fields)
                field = fields{i};
                results(row).(field) = metrics.(field);
            end

            if ~isfield(best_by_mode, mode_key) || is_better(results(row), best_by_mode.(mode_key))
                best_by_mode.(mode_key) = results(row);
            end

            fprintf(['post5L2=%.6f, weighted=%.6f, theta=%.6f, dot=%.6f, ' ...
                     'thetaRMSE=%.4f deg, dotRMSE=%.4f deg/s, maxTheta=%.4f deg, updates=%g\n'], ...
                results(row).combined_l2_5_end, results(row).weighted_l2_5_end, ...
                results(row).theta_l2_5_end, results(row).theta_dot_l2_5_end, ...
                results(row).theta_rmse_5_end, results(row).theta_dot_rmse_5_end, ...
                results(row).max_abs_theta_err_5_end, results(row).kernel_updates);
        catch ME
            results(row).error = ME.message;
            fprintf('FAILED: %s\n', ME.message);
        end
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, {'kernel_mode', 'weighted_l2_5_end', 'combined_l2_5_end'}, ...
    {'ascend', 'ascend', 'ascend'});

if ~exist('Results', 'dir')
    mkdir('Results');
end

save(fullfile('Results', 'long_term_mode2_mode3_retune.mat'), ...
    'results', 'TBL', 'mode_specs', 'best_by_mode');

disp(TBL(:, {'kernel_mode','A_theta','A_theta_dot','Q_theta','Q_theta_dot', ...
    'Gamma_x_theta','Gamma_x_theta_dot','Gamma_r','Gamma_f', ...
    'Gamma_g_theta','Gamma_g_theta_dot','kernel_l2_high','kernel_l2_low', ...
    'weighted_l2_5_end','combined_l2_5_end','theta_l2_5_end', ...
    'theta_dot_l2_5_end','theta_rmse_5_end','theta_dot_rmse_5_end', ...
    'max_abs_theta_err_5_end','kernel_updates','error'}));

end

function C = mode2_candidates()
C = [
     8  9 1.0 1.0 110  70  18  25 110  70 0.020 0.015
    10 12 1.0 1.0 110  70  18  25 110  70 0.020 0.015
    12 14 1.0 2.0 130  80  20  28 130  80 0.018 0.012
    12 16 1.0 3.0 150  90  18  30 150  90 0.018 0.012
    14 18 1.0 5.0 180 100  16  35 180 100 0.016 0.010
    16 20 1.0 5.0 220 120  14  40 220 120 0.016 0.010
    12 16 2.0 5.0 180 100  16  35 180 100 0.016 0.010
    16 22 2.0 8.0 220 130  14  40 220 130 0.014 0.009
    18 24 1.0 8.0 260 150  12  45 260 150 0.014 0.009
    20 26 1.0 10.0 300 170 10  50 300 170 0.012 0.008
    16 22 0.7 6.0 260 140 12  45 260 140 0.014 0.009
    18 24 0.5 8.0 320 180 10  55 320 180 0.012 0.008
    20 28 1.5 12.0 360 200 10  60 360 200 0.012 0.008
    22 30 1.0 15.0 420 240  8  70 420 240 0.010 0.006
];
end

function C = mode3_candidates()
C = [
     8  9 1.0 1.0  70  40  26 17  70  40 0.020 0.015
    10 12 1.0 1.0 120  70  24 22 120  70 0.020 0.015
    12 14 1.0 2.0 180 100  22 28 180 100 0.018 0.012
    12 16 1.0 3.0 260 140  20 35 260 140 0.018 0.012
    14 18 1.0 5.0 340 180  18 45 340 180 0.016 0.010
    16 20 1.0 5.0 450 240  16 55 450 240 0.016 0.010
    18 24 1.0 8.0 600 320  14 70 600 320 0.014 0.009
    20 26 1.0 10.0 750 420 12 85 750 420 0.012 0.008
    18 24 2.0 10.0 600 320 14 70 600 320 0.014 0.009
    20 28 2.0 15.0 750 420 12 85 750 420 0.012 0.008
    18 24 0.7 8.0 750 380 12 80 750 380 0.012 0.008
    22 30 1.0 15.0 900 500 10 100 900 500 0.010 0.006
    24 32 1.0 20.0 1100 620 8 120 1100 620 0.010 0.006
    24 34 0.5 20.0 1300 700 8 140 1300 700 0.009 0.005
];
end

function overrides = make_overrides(kernel_mode, cand)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', kernel_mode, ...
    'plant_schedule', 'aggressive', ...
    'Tmax', 10, ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'A_tran', diag(-cand(1:2)), ...
    'Q', diag(cand(3:4)), ...
    'Gamma_x', diag(cand(5:6)), ...
    'Gamma_r', cand(7), ...
    'Gamma_f', cand(8), ...
    'Gamma_g', diag(cand(9:10)), ...
    'kernel_l2_high', cand(11), ...
    'kernel_l2_low', cand(12), ...
    'kernel_rep_point', 'endpoint', ...
    'max_depth', 7);
end

function metrics = run_one_candidate(overrides)
setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:,1:2) - X(:,3:4);
idx_tail = T >= 5;
if nnz(idx_tail) < 2
    error('Not enough samples after 5 seconds.');
end

t_tail = T(idx_tail);
tau = (t_tail - 5) ./ max(t_tail(end) - 5, eps);
w = 1 + 2*tau;
tail_err = e(idx_tail, :);
dot_err = tail_err(:, 2);

metrics = struct();
metrics.theta_l2_5_end = sqrt(trapz(t_tail, tail_err(:,1).^2));
metrics.theta_dot_l2_5_end = sqrt(trapz(t_tail, dot_err.^2));
metrics.combined_l2_5_end = sqrt(trapz(t_tail, sum(tail_err.^2, 2)));
metrics.weighted_l2_5_end = sqrt(trapz(t_tail, w .* sum(tail_err.^2, 2)));
metrics.theta_rmse_5_end = sqrt(mean(rad2deg(tail_err(:,1)).^2));
metrics.theta_dot_rmse_5_end = sqrt(mean(rad2deg(dot_err).^2));
metrics.max_abs_theta_err_5_end = max(abs(rad2deg(tail_err(:,1))));
if exist('SIM_METRICS', 'var') == 1 && isfield(SIM_METRICS, 'control_rms')
    metrics.control_rms = SIM_METRICS.control_rms;
    metrics.control_max_abs = SIM_METRICS.control_max_abs;
else
    metrics.control_rms = NaN;
    metrics.control_max_abs = NaN;
end
if exist('SIM_METRICS', 'var') == 1 && isfield(SIM_METRICS, 'kernel_updates')
    metrics.kernel_updates = SIM_METRICS.kernel_updates;
else
    metrics.kernel_updates = numel(center_times);
end
end

function tf = is_better(a, b)
tf = a.weighted_l2_5_end < b.weighted_l2_5_end;
if ~tf && abs(a.weighted_l2_5_end - b.weighted_l2_5_end) < 1e-12
    tf = a.combined_l2_5_end < b.combined_l2_5_end;
end
end

function result = empty_result()
result = struct( ...
    'kernel_mode', NaN, ...
    'A_theta', NaN, ...
    'A_theta_dot', NaN, ...
    'Q_theta', NaN, ...
    'Q_theta_dot', NaN, ...
    'Gamma_x_theta', NaN, ...
    'Gamma_x_theta_dot', NaN, ...
    'Gamma_r', NaN, ...
    'Gamma_f', NaN, ...
    'Gamma_g_theta', NaN, ...
    'Gamma_g_theta_dot', NaN, ...
    'kernel_l2_high', NaN, ...
    'kernel_l2_low', NaN, ...
    'theta_l2_5_end', NaN, ...
    'theta_dot_l2_5_end', NaN, ...
    'combined_l2_5_end', NaN, ...
    'weighted_l2_5_end', NaN, ...
    'theta_rmse_5_end', NaN, ...
    'theta_dot_rmse_5_end', NaN, ...
    'max_abs_theta_err_5_end', NaN, ...
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
