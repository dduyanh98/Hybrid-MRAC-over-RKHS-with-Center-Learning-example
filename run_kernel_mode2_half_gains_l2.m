function TBL = run_kernel_mode2_half_gains_l2()
%RUN_KERNEL_MODE2_HALF_GAINS_L2 Run architecture 3 / kernel mode 2 at half gains.

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
mode2 = p0.kernel_mode2_params;

overrides = struct( ...
    'architecture_mode', 3, ...
    'kernel_mode', 2, ...
    'Gamma_x', 0.5 * mode2.Gamma_x, ...
    'Gamma_r', 0.5 * mode2.Gamma_r, ...
    'Gamma_f', 0.5 * mode2.Gamma_f, ...
    'Gamma_g', 0.5 * mode2.Gamma_g);

setenv('PENDULUM_OVERRIDES', jsonencode(overrides));
evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
e = X(:, 1:2) - X(:, 3:4);

controller = "Kernel mode 2 with half current gains";
kernel_mode = 2;
gain_scale = 0.5;
Gamma_x_theta = overrides.Gamma_x(1, 1);
Gamma_x_theta_dot = overrides.Gamma_x(2, 2);
Gamma_r = overrides.Gamma_r;
Gamma_f = overrides.Gamma_f;
Gamma_g_theta = overrides.Gamma_g(1, 1);
Gamma_g_theta_dot = overrides.Gamma_g(2, 2);
L2_theta = sqrt(trapz(T, e(:, 1).^2));
L2_theta_dot = sqrt(trapz(T, e(:, 2).^2));
L2_x_minus_xref = sqrt(trapz(T, sum(e.^2, 2)));
kernel_updates = SIM_METRICS.kernel_updates;

TBL = table(controller, kernel_mode, gain_scale, ...
    Gamma_x_theta, Gamma_x_theta_dot, Gamma_r, Gamma_f, ...
    Gamma_g_theta, Gamma_g_theta_dot, L2_theta, L2_theta_dot, ...
    L2_x_minus_xref, kernel_updates);

disp(TBL);
save(fullfile(out_dir, 'kernel_mode2_half_gains_l2.mat'), 'TBL');
writetable(TBL, fullfile(out_dir, 'kernel_mode2_half_gains_l2.csv'));

fprintf('Saved table to %s\n', out_dir);
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
