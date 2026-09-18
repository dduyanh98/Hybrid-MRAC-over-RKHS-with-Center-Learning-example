function TBL = run_architecture_kernel_l2_table()
%RUN_ARCHITECTURE_KERNEL_L2_TABLE Compare tracking-error L2 norms.
% Rows:
%   1) Architecture 2
%   2) Architecture 3, kernel center mode 2
%   3) Architecture 3, kernel center mode 3

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

cases = [
    struct('label', "Architecture 2", ...
           'architecture_mode', 2, ...
           'kernel_mode', NaN)
    struct('label', "Architecture 3, kernel mode 2", ...
           'architecture_mode', 3, ...
           'kernel_mode', 2)
    struct('label', "Architecture 3, kernel mode 3", ...
           'architecture_mode', 3, ...
           'kernel_mode', 3)
];

controller = strings(numel(cases), 1);
architecture_mode = zeros(numel(cases), 1);
kernel_mode = nan(numel(cases), 1);
L2_theta = nan(numel(cases), 1);
L2_theta_dot = nan(numel(cases), 1);

for i = 1:numel(cases)
    fprintf('Running %s...\n', cases(i).label);
    metrics = run_one_case(cases(i));

    controller(i) = cases(i).label;
    architecture_mode(i) = cases(i).architecture_mode;
    kernel_mode(i) = cases(i).kernel_mode;
    L2_theta(i) = metrics.L2_theta;
    L2_theta_dot(i) = metrics.L2_theta_dot;

    fprintf('  L2 theta = %.6g, L2 theta_dot = %.6g\n', ...
        L2_theta(i), L2_theta_dot(i));
end

TBL = table(controller, architecture_mode, kernel_mode, L2_theta, L2_theta_dot);

disp(TBL);
save(fullfile(out_dir, 'architecture_kernel_l2_table.mat'), 'TBL');
writetable(TBL, fullfile(out_dir, 'architecture_kernel_l2_table.csv'));

fprintf('Saved table to %s\n', out_dir);
end

function metrics = run_one_case(test_case)
overrides = struct('architecture_mode', test_case.architecture_mode);
if ~isnan(test_case.kernel_mode)
    overrides.kernel_mode = test_case.kernel_mode;
end

setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:, 1:2) - X(:, 3:4);

metrics = struct();
metrics.L2_theta = sqrt(trapz(T, e(:, 1).^2));
metrics.L2_theta_dot = sqrt(trapz(T, e(:, 2).^2));
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
