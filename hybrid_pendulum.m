function dy = hybrid_pendulum(t, y, p, Xi, alpha)
% HYBRID_PENDULUM — supports architectures 1–3.
% State layout:
%   1:2   x          (plant states)
%   3:4   x_ref
%   5:6   K_hat_x
%   7     K_hat_r
%   8     I_ref
%   9:10  e_tran
%   11    Ieps      = ∫ W(ε) dt
%   12    Itran     = ∫ W(e_tran) dt
%   13    Ex_f_hat  (adaptive-only in mode 1)
%   14:15 K_hat_g
%   p.idx_theta_rkhs  RKHS weights
%   p.idx_S_eps       jump sum for eps'P eps
%   p.idx_S_tran      jump sum for e_tran'P_tran e_tran
%   p.idx_kernel_l2   integral of ||x - x_ref||^2 for kernel events

% ------------------------------------------------------------------------
% unpack states
% ------------------------------------------------------------------------
x        = y(1:2);
x_ref    = y(3:4);
K_hat_x  = y(5:6);
K_hat_r  = y(7);
I_ref    = y(8);
e_tran   = y(9:10);
Ieps     = y(11); %#ok<NASGU>
Itran    = y(12); %#ok<NASGU>
Theta_hat = y(13);
K_hat_g        = y(14:15);
Theta_rkhs     = y(p.idx_theta_rkhs);

% ------------------------------------------------------------------------
% parameters
% ------------------------------------------------------------------------
m = p.m;     l = p.l;     k = p.k;     g = p.g;
m_est = p.m_est; l_est = p.l_est; k_est = p.k_est;
kp = p.kp; kd = p.kd; ki = p.ki;
Gamma_x = p.Gamma_x; Gamma_r = p.Gamma_r; Gamma_f = p.Gamma_f; Gamma_g = p.Gamma_g;
P = p.P; Q = p.Q; A_tran = p.A_tran;

% we also have V_fun, W_fun from hybrid_params; we will use W_fun
if ~isfield(p,'W_fun') || isempty(p.W_fun)
    % fallback: W(z) = zᵀQz
    p.W_fun = @(z) z(:)'*Q*z(:);
end

% ------------------------------------------------------------------------
% reference model
% ------------------------------------------------------------------------
yd   = p.desired_output(t);
yd_d = p.desired_output_dot(t);

A_ref = [0 1; -kp -kd];
B_ref = [0; 1];

r     = kp*yd - ki*I_ref;
x_ref_dot = A_ref*x_ref + B_ref*r;

% ------------------------------------------------------------------------
% main error channels
% ------------------------------------------------------------------------
e   = x - x_ref;        % tracking error
eps = e - e_tran;       % hybrid-layer epsilon

% ------------------------------------------------------------------------
% true regressor (plant uncertainty)
% ------------------------------------------------------------------------
phi_true = regressor_vector(x);
Theta_true = -m*g;
Ex_f = Theta_true * phi_true;

% ------------------------------------------------------------------------
% estimated Ex_f_hat (adaptive-only or kernel)
% ------------------------------------------------------------------------
use_kernel = (p.architecture_mode >= 2) && ~isempty(Xi) && ~isempty(alpha);

if use_kernel
    % RKHS approximation: Theta_rkhs' * Phi_rkhs(x).
    [~, K_aprx] = kernel_prediction(x', Xi, zeros(size(alpha,1),1), p.sigma, 1);
    Phi_rkhs = K_aprx(:);
    n_rkhs = min(numel(Theta_rkhs), numel(Phi_rkhs));
    Ex_f_hat = Theta_rkhs(1:n_rkhs).' * Phi_rkhs(1:n_rkhs);
else
    % adaptive-only scalar (mode 1)
    Phi_rkhs = zeros(4,1);
    Ex_f_hat = Theta_hat * phi_true;
end

PB = P * [0; 1/(m*l)];
if p.architecture_mode == 3
    adapt_error = eps;
else
    adapt_error = e;
end
sigma_ad = adapt_error.' * PB;

Theta_hat_dot = Gamma_f * (phi_true * sigma_ad);

% ------------------------------------------------------------------------
% control law
% ------------------------------------------------------------------------
theta     = x(1);
theta_dot = x(2);

u_baseline = m_est*g*sin(theta) + k_est*l_est*theta_dot ...
            - m_est*l_est*( kp*(theta - yd) + kd*(theta_dot - yd_d) );

u_adaptive = K_hat_x.' * x + K_hat_r * r - Ex_f_hat;
if p.architecture_mode == 3
    u_adaptive = u_adaptive + K_hat_g.' * e;
end

u = u_baseline + u_adaptive;

% ------------------------------------------------------------------------
% plant dynamics
% ------------------------------------------------------------------------
A = [0 1;
     0 -k/m];
B = [0;
     1/(m*l)];
x_dot = A*x + B*(u + Ex_f);

% ------------------------------------------------------------------------
% adaptive laws
% ------------------------------------------------------------------------
K_hat_x_dot = -Gamma_x * (x * sigma_ad);
K_hat_r_dot = -Gamma_r * (r * sigma_ad);
K_hat_g_dot = zeros(2,1);
if p.architecture_mode == 3
    K_hat_g_dot = -Gamma_g * (e * sigma_ad);
end

Theta_rkhs_dot = zeros(numel(Theta_rkhs),1);
if use_kernel
    n_rkhs = min(numel(Theta_rkhs), numel(Phi_rkhs));
    Theta_rkhs_dot(1:n_rkhs) = Gamma_f * ...
        (alpha(1:n_rkhs,1:n_rkhs) * Phi_rkhs(1:n_rkhs)) * sigma_ad;
end

% ------------------------------------------------------------------------
% auxiliary integrals (Eq. 46, 48)
%   Ieps  = ∫ W(ε)
%   Itran = ∫ W(e_tran)
% ------------------------------------------------------------------------
I_ref_dot  = x_ref(1) - yd;      % Soonyong-style integral (keep as you had)

% (46)  W(ε)
Ieps_dot  = p.W_fun(eps);

% (48)  W(e_tran)
Itran_dot = p.W_fun(e_tran);

% ------------------------------------------------------------------------
% transient channel dynamics
% ------------------------------------------------------------------------
e_tran_dot = A_tran * e_tran;

% Eq. (46), (48) use sums of Lyapunov jumps. Main_sim updates these states
% at hybrid events; they stay constant during continuous flow.
S_eps_dot = 0;
S_tran_dot = 0;
kernel_l2_dot = 0;
if p.architecture_mode == 3 && isfield(p,'idx_kernel_l2')
    kernel_l2_dot = e.' * e;
end

% ------------------------------------------------------------------------
% pack derivative vector
% ------------------------------------------------------------------------
dy = zeros(size(y));

dy(1:2)   = x_dot;
dy(3:4)   = x_ref_dot;
dy(5:6)   = K_hat_x_dot;
dy(7)     = K_hat_r_dot;
dy(8)     = I_ref_dot;

dy(9:10)  = e_tran_dot;
dy(11)    = Ieps_dot;
dy(12)    = Itran_dot;

dy(13)    = Theta_hat_dot;
dy(14:15) = K_hat_g_dot;
dy(p.idx_theta_rkhs) = Theta_rkhs_dot;
dy(p.idx_S_eps)      = S_eps_dot;
dy(p.idx_S_tran)     = S_tran_dot;
if isfield(p,'idx_kernel_l2')
    dy(p.idx_kernel_l2) = kernel_l2_dot;
end

% ------------------------------------------------------------------------
% debug print: norms and Ieps_dot every ~0.1s
% ------------------------------------------------------------------------
if isfield(p,'debug') && p.debug && mod(t, 0.1) < 1e-4
    e_norm      = norm(e);
    e_tran_norm = norm(e_tran);
    eps_norm    = norm(eps);
    fprintf('t = %.3f: [‖e‖ %.4f  ‖e_tran‖ %.4f  ‖eps‖ %.4f  Ieps_dot %.4g]\n', ...
        t, e_norm, e_tran_norm, eps_norm, Ieps_dot);
end

end
