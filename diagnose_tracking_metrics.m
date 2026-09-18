function metrics = diagnose_tracking_metrics()
%DIAGNOSE_TRACKING_METRICS Compare model-reference and desired-output errors.

addpath(genpath(pwd));

old_no_plots = getenv('PENDULUM_NO_PLOTS');
old_overrides = getenv('PENDULUM_OVERRIDES');
old_keep_workspace = getenv('PENDULUM_KEEP_WORKSPACE');
cleanup = onCleanup(@() restore_env(old_no_plots, old_overrides, old_keep_workspace));

setenv('PENDULUM_NO_PLOTS', '1');
setenv('PENDULUM_KEEP_WORKSPACE', '1');
setenv('PENDULUM_OVERRIDES', '');

evalc('run(''Main_sim.m'')');

T = LOG.T;
X = LOG.X;
x = X(:,1:2);
xref = X(:,3:4);
e_ref = x - xref;

yd = arrayfun(p.desired_output, T);
yd_dot = arrayfun(p.desired_output_dot, T);
e_des = [x(:,1) - yd(:), x(:,2) - yd_dot(:)];
e_ref_des = [xref(:,1) - yd(:), xref(:,2) - yd_dot(:)];

idx_late = T >= 5 & T <= 10;

metrics = struct();
metrics.model_ref_l2_5_10_raw = sqrt(trapz(T(idx_late), sum(e_ref(idx_late,:).^2, 2)));
metrics.desired_l2_5_10_raw = sqrt(trapz(T(idx_late), sum(e_des(idx_late,:).^2, 2)));
metrics.ref_model_to_desired_l2_5_10_raw = sqrt(trapz(T(idx_late), sum(e_ref_des(idx_late,:).^2, 2)));

metrics.model_ref_theta_rmse_deg = sqrt(mean(rad2deg(e_ref(idx_late,1)).^2));
metrics.model_ref_theta_dot_rmse_deg_s = sqrt(mean(rad2deg(e_ref(idx_late,2)).^2));
metrics.desired_theta_rmse_deg = sqrt(mean(rad2deg(e_des(idx_late,1)).^2));
metrics.desired_theta_dot_rmse_deg_s = sqrt(mean(rad2deg(e_des(idx_late,2)).^2));
metrics.ref_model_to_desired_theta_rmse_deg = sqrt(mean(rad2deg(e_ref_des(idx_late,1)).^2));
metrics.ref_model_to_desired_theta_dot_rmse_deg_s = sqrt(mean(rad2deg(e_ref_des(idx_late,2)).^2));

metrics.model_ref_theta_max_deg = max(abs(rad2deg(e_ref(idx_late,1))));
metrics.desired_theta_max_deg = max(abs(rad2deg(e_des(idx_late,1))));
metrics.ref_model_to_desired_theta_max_deg = max(abs(rad2deg(e_ref_des(idx_late,1))));

metrics.P_e_l2_5_10 = sqrt(trapz(T(idx_late), ...
    arrayfun(@(i) e_ref(i,:) * p.P * e_ref(i,:)', find(idx_late)).^2));

disp(metrics);

save(fullfile('Results', 'tracking_metric_diagnosis.mat'), ...
    'metrics', 'T', 'x', 'xref', 'yd', 'yd_dot', 'e_ref', 'e_des', 'e_ref_des');
end

function restore_env(old_no_plots, old_overrides, old_keep_workspace)
setenv('PENDULUM_NO_PLOTS', old_no_plots);
setenv('PENDULUM_OVERRIDES', old_overrides);
setenv('PENDULUM_KEEP_WORKSPACE', old_keep_workspace);
end
