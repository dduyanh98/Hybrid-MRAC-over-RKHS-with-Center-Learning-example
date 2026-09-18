function TBL = run_current_setting_gain_tune()
%RUN_CURRENT_SETTING_GAIN_TUNE Tune adaptive gains for current Main_sim setup.
%
% Fixed setup:
%   - architecture_mode = 3
%   - kernel_mode = 3
%   - plant_schedule = current
%   - fixed/constant hybrid budgets
%   - endpoint representative
%   - kernel window high/low = 0.02/0.015
%
% Ranking metric: L2 tracking norm over 5-10 seconds.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

gain_table = [
     20   20   10   20
     50   50   40   80
    100  150   50  100
    200  220   75  200
    250  260  100  250
    300  320  125  300
    300  350  150  300
    400  350  150  400
    500  340  175  500
    600  450  225  600
];

results = repmat(empty_result(), size(gain_table, 1), 1);
best_run = struct();

for i = 1:size(gain_table, 1)
    gains = gain_table(i, :);
    fprintf('\n===== CURRENT SETTING GAIN TEST %d/%d: gx=%g gr=%g gf=%g gg=%g =====\n', ...
        i, size(gain_table, 1), gains(1), gains(2), gains(3), gains(4));

    results(i).Gamma_x = gains(1);
    results(i).Gamma_r = gains(2);
    results(i).Gamma_f = gains(3);
    results(i).Gamma_g = gains(4);

    try
        run_data = run_one_candidate(make_overrides(gains));
        fields = fieldnames(run_data.metrics);
        for f = 1:numel(fields)
            field = fields{f};
            results(i).(field) = run_data.metrics.(field);
        end

        if ~isfield(best_run, 'metrics') || ...
                run_data.metrics.l2_5_10 < best_run.metrics.l2_5_10
            best_run = run_data;
            best_run.gains = gains;
        end

        fprintf(['L2_5_10=%.6f, theta=%.4f deg, thetaDot=%.4f deg/s, ' ...
                 'maxTheta=%.4f deg, updates=%g, uRMS=%.4f, uMax=%.4f\n'], ...
            results(i).l2_5_10, results(i).theta_rmse_5_10, ...
            results(i).theta_dot_rmse_5_10, results(i).max_abs_theta_err_5_10, ...
            results(i).kernel_updates, results(i).control_rms, results(i).control_max_abs);
    catch ME
        results(i).error = ME.message;
        fprintf('FAILED: %s\n', ME.message);
    end
end

TBL = struct2table(results);
TBL = sortrows(TBL, 'l2_5_10', 'ascend');

if ~exist('Results', 'dir')
    mkdir('Results');
end

if isfield(best_run, 'metrics')
    plot_best_run(best_run);
end

save(fullfile('Results', 'current_setting_gain_tune.mat'), ...
    'results', 'TBL', 'gain_table', 'best_run');

disp(TBL);
end

function overrides = make_overrides(gains)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', 3, ...
    'plant_schedule', 'current', ...
    'hybrid_budget_mode', 'fixed', ...
    'fixed_budget_law', 'constant', ...
    'kernel_rep_point', 'endpoint', ...
    'kernel_l2_high', 0.02, ...
    'kernel_l2_low', 0.015, ...
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

metrics = struct();
metrics.l2_5_10 = sqrt(trapz(T(idx_late), sum(e(idx_late,:).^2, 2)));
metrics.l2_full = SIM_METRICS.tracking_l2;
metrics.theta_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,1)).^2));
metrics.theta_dot_rmse_5_10 = sqrt(mean(rad2deg(e(idx_late,2)).^2));
metrics.max_abs_theta_err_5_10 = max(abs(rad2deg(e(idx_late,1))));
metrics.control_rms = SIM_METRICS.control_rms;
metrics.control_max_abs = SIM_METRICS.control_max_abs;
metrics.kernel_updates = SIM_METRICS.kernel_updates;

run_data = struct('T', T, 'e', e, 'metrics', metrics);
end

function plot_best_run(best_run)
figure('Color', 'w', 'Visible', 'off');
subplot(2,1,1);
plot(best_run.T, rad2deg(best_run.e(:,1)), 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('e_\theta [deg]');
grid on;
title(sprintf('Best Current Gain Error: G=[%g %g %g %g], L2_{5-10}=%.4g', ...
    best_run.gains(1), best_run.gains(2), best_run.gains(3), best_run.gains(4), ...
    best_run.metrics.l2_5_10));

subplot(2,1,2);
plot(best_run.T, rad2deg(best_run.e(:,2)), 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('e_{\dot{\theta}} [deg/s]');
grid on;

saveas(gcf, fullfile('Results', 'current_setting_gain_tune_best_error.png'));
savefig(gcf, fullfile('Results', 'current_setting_gain_tune_best_error.fig'));
end

function result = empty_result()
result = struct( ...
    'Gamma_x', NaN, ...
    'Gamma_r', NaN, ...
    'Gamma_f', NaN, ...
    'Gamma_g', NaN, ...
    'l2_5_10', NaN, ...
    'l2_full', NaN, ...
    'theta_rmse_5_10', NaN, ...
    'theta_dot_rmse_5_10', NaN, ...
    'max_abs_theta_err_5_10', NaN, ...
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
