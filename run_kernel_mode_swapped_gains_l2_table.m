function TBL = run_kernel_mode_swapped_gains_l2_table()
%RUN_KERNEL_MODE_SWAPPED_GAINS_L2_TABLE Compare kernel modes with swapped gains.
% Runs architecture 3 twice:
%   1) kernel mode 2 using kernel mode 3 adaptive gains
%   2) kernel mode 3 using kernel mode 2 adaptive gains

addpath(genpath(pwd));

out_dir = fullfile(pwd, 'Results');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');

p0 = hybrid_params();
mode2_gains = gain_fields(p0.kernel_mode2_params);
mode3_gains = gain_fields(p0.kernel_mode3_params);

cases = [
    struct('controller', "Kernel mode 2 with mode 3 gains", ...
           'kernel_mode', 2, ...
           'gains', mode3_gains)
    struct('controller', "Kernel mode 3 with mode 2 gains", ...
           'kernel_mode', 3, ...
           'gains', mode2_gains)
];

controller = strings(numel(cases), 1);
kernel_mode = zeros(numel(cases), 1);
gain_source = strings(numel(cases), 1);
L2_theta = nan(numel(cases), 1);
L2_theta_dot = nan(numel(cases), 1);
L2_x_minus_xref = nan(numel(cases), 1);
kernel_updates = nan(numel(cases), 1);

for i = 1:numel(cases)
    fprintf('Running %s...\n', cases(i).controller);
    metrics = run_one_case(cases(i));

    controller(i) = cases(i).controller;
    kernel_mode(i) = cases(i).kernel_mode;
    gain_source(i) = sprintf('kernel mode %d gains', 5 - cases(i).kernel_mode);
    L2_theta(i) = metrics.L2_theta;
    L2_theta_dot(i) = metrics.L2_theta_dot;
    L2_x_minus_xref(i) = metrics.L2_x_minus_xref;
    kernel_updates(i) = metrics.kernel_updates;

    fprintf('  L2 theta = %.6g, L2 theta_dot = %.6g, L2 x-xref = %.6g\n', ...
        L2_theta(i), L2_theta_dot(i), L2_x_minus_xref(i));
end

TBL = table(controller, kernel_mode, gain_source, L2_theta, ...
    L2_theta_dot, L2_x_minus_xref, kernel_updates);

disp(TBL);
save(fullfile(out_dir, 'kernel_mode_swapped_gains_l2_table.mat'), 'TBL');
writetable(TBL, fullfile(out_dir, 'kernel_mode_swapped_gains_l2_table.csv'));

fprintf('Saved table to %s\n', out_dir);
end

function gains = gain_fields(params)
gains = struct( ...
    'Gamma_x', params.Gamma_x, ...
    'Gamma_r', params.Gamma_r, ...
    'Gamma_f', params.Gamma_f, ...
    'Gamma_g', params.Gamma_g);
end

function metrics = run_one_case(test_case)
overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', test_case.kernel_mode, ...
    'Gamma_x', test_case.gains.Gamma_x, ...
    'Gamma_r', test_case.gains.Gamma_r, ...
    'Gamma_f', test_case.gains.Gamma_f, ...
    'Gamma_g', test_case.gains.Gamma_g);

setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:, 1:2) - X(:, 3:4);

metrics = struct();
metrics.L2_theta = sqrt(trapz(T, e(:, 1).^2));
metrics.L2_theta_dot = sqrt(trapz(T, e(:, 2).^2));
metrics.L2_x_minus_xref = sqrt(trapz(T, sum(e.^2, 2)));
metrics.kernel_updates = SIM_METRICS.kernel_updates;
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
