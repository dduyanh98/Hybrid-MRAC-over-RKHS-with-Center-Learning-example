function p = hybrid_params()
% HYBRID_PARAMS — Parameter struct for hybrid pendulum (arch 1–3)
% Clean, reproducible, no persistent or stale fields.

%% --- Physical parameters (Initial plant) ----------------------------------
p.m = 1.3;              % mass [kg]
p.l = 0.8;              % length [m]
p.k = 0.5;              % viscous damping
p.g = 9.81;             % gravity

%% --- Estimated parameters (baseline model) -----------------------------
p.m_est = 1.0;
p.l_est = 1.0;
p.k_est = 0.1;

%% --- Reference model gains ---------------------------------------------
wn   = sqrt(p.g / p.l_est);
zeta = 0.75;

p.kp = wn^2;
p.kd = 2*zeta*wn;
p.ki = 0.1;

%% --- Adaptive gains ----------------------------------------------------
p.Gamma_x = 5 * eye(2);
p.Gamma_r = 5;
p.Gamma_f = 4;
p.Gamma_g = 8 * eye(2);

%% --- Lyapunov matrices -------------------------------------------------
p.Q = eye(2);
lambda_min_Q = min(eig(p.Q));
p.W_fun = @(z) lambda_min_Q * (z(:)' * z(:));

%% --- Transient channel -------------------------------------------------
p.A_tran = diag([-8, -9]);
p.Q_tran = eye(2);      
p.P = lyap(p.A_tran', p.Q);
p.P_tran = lyap(p.A_tran', p.Q_tran);

%% Eq. (52)-(53) reset shaping maps. With these default choices,
% P^(-1/2) h_ref = epsilon and P_tran^(-1/2) h_tran = e_tran.
p.h_ref = @(t, eps, P) sqrtm(P) * eps;
p.h_tran = @(t, eps, P_tran) sqrtm(P_tran) * eps;

%% --- Dictionary parameters ---------------------------------------------
p.r0        = 5.0;
p.gamma     = -0.1;
p.rho       = 0.25;
p.max_depth = 4;


%% --- Kernel layer (arch >= 2) -----------------------------------------
p.sigma  = 0.75;
p.r_min  = 0.10;
p.r_grid = 0.10;
p.kernel_l2_high = 0.01;
p.kernel_l2_low  = 0.005;
p.kernel_rep_point = 'endpoint'; % 'endpoint' or 'max_tracking_error'
p.kernel_mode1_min_grid_n = 2;
p.kernel_mode1_max_grid_n = 5;
p.kernel_mode2_min_grid_n = 2;
p.kernel_mode2_max_grid_n = 6;
p.show_center_change_lines = false;
p.kernel_min_centers = p.kernel_mode2_min_grid_n^2;
p.kernel_max_centers = 40;
p.idx_theta_rkhs = 16:(15+p.kernel_max_centers);
p.idx_S_eps = 16+p.kernel_max_centers;
p.idx_S_tran = 17+p.kernel_max_centers;
p.idx_kernel_l2 = 18+p.kernel_max_centers;

%% --- Kernel-mode-specific tuned parameters -----------------------------
% These are applied in Main_sim after p.kernel_mode is selected.
p.kernel_mode2_params = struct( ...
    'Gamma_x', diag([110, 70]), ...
    'Gamma_r', 18, ...
    'Gamma_f', 25, ...
    'Gamma_g', diag([110, 70]), ...
    'Q', eye(2), ...
    'kernel_l2_high', 0.02, ...
    'kernel_l2_low', 0.015, ...
    'max_depth', 7);

p.kernel_mode3_params = struct( ...
    'Gamma_x', diag([70, 40]), ...
    'Gamma_r', 26, ...
    'Gamma_f', 17, ...
    'Gamma_g', diag([70, 40]), ...
    'Q', diag([1, 1]), ...
    'kernel_l2_high', 0.02, ...
    'kernel_l2_low', 0.015, ...
    'kernel_rep_point', 'endpoint', ...
    'max_depth', 7);

%% --- Epoch & hybrid thresholds -----------------------------------------
p.reset_series_time = 0.1;
p.plant_schedule = 'aggressive';     % 'current' or 'aggressive'

%% --- Hybrid Lyapunov bookkeeping (RESET EVERY RUN) ---------------------
p.ref_budget_cum  = 0;
p.tran_budget_cum = 0;
p.event_tol = 1e-10;
p.alpha_hybrid_series = 1.001;
p.hybrid_budget_mode = 'fixed';  % 'faithful' or 'fixed'
p.fixed_budget_law = 'constant'; % 'constant', 'floor'/'floor_actual', or 'positive'/'positive_actual'
p.positive_jump_budget_only = true;
p.accumulate_flow_budget = false;
p.fixed_integral_jump_correction = false;
p.DeltaV_step_ref  = 1e-3;
p.DeltaW_step_tran = 1e-3;
p.min_jump_budget_eps = 0;
p.min_jump_budget_tran = 0;

p.s_ref  = 0;
p.s_tran = 0;

p.Ieps_anchor  = 0;
p.Itran_anchor = 0;
p.S_eps_ref_anchor = 0;
p.S_eps_tran_anchor = 0;
p.S_tran_anchor = 0;
p.budget_tol = 1e-12;
p.event_time_eps = 1e-4;
p.t_last_ref_jump = 0;
p.t_last_tran_jump = 0;

p.ref_jump_times    = [];
p.tran_jump_times   = [];
p.budget_jump_times = [];
p.DeltaV_hist       = [];
p.DeltaW_hist       = [];
p.Ieps_anchor_hist  = [];
p.Itran_anchor_hist = [];

%% --- Epoch-time bookkeeping -------------------------------------------
p.t_last_epoch = 0;

end
