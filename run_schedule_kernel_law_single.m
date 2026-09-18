function result = run_schedule_kernel_law_single(schedule, case_name)
%RUN_SCHEDULE_KERNEL_LAW_SINGLE Run one schedule/kernel-law candidate.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

cases = schedule_kernel_law_case_defs();
case_idx = find(strcmp({cases.name}, case_name), 1);
if isempty(case_idx)
    error('Unknown case_name: %s', case_name);
end

test_case = cases(case_idx);
overrides = make_overrides(schedule, test_case.params);

fprintf('Running %s / %s\n', schedule, case_name);
result = empty_result();
result.schedule = schedule;
result.case_name = case_name;
result.kernel_mode = test_case.params.kernel_mode;
result.kernel_l2_high = test_case.params.kernel_l2_high;
result.kernel_l2_low = test_case.params.kernel_l2_low;

try
    run_data = run_one_candidate(overrides);
    fields = fieldnames(run_data.metrics);
    for i = 1:numel(fields)
        field = fields{i};
        result.(field) = run_data.metrics.(field);
    end
    result.T = run_data.T;
    result.e = run_data.e;
    result.error = '';
catch ME
    result.error = ME.message;
    result.T = [];
    result.e = [];
end

out_dir = fullfile('Results', 'schedule_kernel_law_cases');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save(fullfile(out_dir, [safe_name(schedule) '__' safe_name(case_name) '.mat']), 'result');

if isempty(result.error)
    fprintf('Done: L2_5_10=%.6f, theta=%.4f deg, thetaDot=%.4f deg/s\n', ...
        result.l2_5_10, result.theta_rmse_5_10, result.theta_dot_rmse_5_10);
else
    fprintf('Failed: %s\n', result.error);
end
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

run_data = struct('T', T, 'e', e, 'metrics', metrics);
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
    'T', [], ...
    'e', [], ...
    'error', '');
end

function name = safe_name(name)
name = regexprep(name, '[^A-Za-z0-9_]+', '_');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
