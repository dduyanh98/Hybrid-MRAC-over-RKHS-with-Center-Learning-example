function [y_plus, p] = apply_jump(which_jump, t, y, p)
% APPLY_JUMP - Jump map for epoch, reference, and transient resets.
%   which_jump = 1 : epoch jump
%   which_jump = 2 : reference jump
%   which_jump = 3 : transient jump

x      = y(1:2);
x_ref  = y(3:4);
e_tran = y(9:10);
Ieps   = y(11);
Itran  = y(12);
S_eps  = y(p.idx_S_eps);
S_tran = y(p.idx_S_tran);

e   = x - x_ref;
eps = e - e_tran;

y_plus = y;

W_eps = eps' * p.P * eps;
W_et  = e_tran' * p.P_tran * e_tran;

if ~isfield(p,'W_fun') || isempty(p.W_fun)
    p.W_fun = @(z) z(:)' * p.Q * z(:);
end

if isfield(p,'alpha_hybrid_series') && ~isempty(p.alpha_hybrid_series)
    alpha_series = p.alpha_hybrid_series;
else
    alpha_series = 1.001;
end

switch which_jump
    case 1
        p = update_plant_at_epoch(p);
        p.t_last_epoch = t;
        return;

    case 2
        if W_eps <= 0
            return;
        end

        if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
            sigma = min(p.DeltaV_step_ref, W_eps);
            alpha_fixed = 1 - sqrt(max(0, (W_eps - sigma)) / W_eps);
            x_ref_plus = x - e_tran - alpha_fixed * eps;
            y_plus(3:4) = x_ref_plus;

            if isfield(p,'fixed_integral_jump_correction') && p.fixed_integral_jump_correction
                e_plus = x - x_ref_plus;
                eps_plus = e_plus - e_tran;
                Ieps_jump_plus = eps_plus' * p.P * eps_plus;
                y_plus(11) = max(0, Ieps + (Ieps_jump_plus - W_eps));
            end

            p.Ieps_anchor = y_plus(11);
            p.ref_jump_times(end+1,1) = t;
            p.Ieps_anchor_hist(end+1,1) = y_plus(11);
            return;
        end

        p.s_ref = find_s(p.s_ref, W_eps, alpha_series);
        sigma = min(series_element(p.s_ref, alpha_series), W_eps);

        P_half = sqrtm(p.P);
        h_ref = p.h_ref(t, eps, p.P);
        h_ref_sq = h_ref' * h_ref;
        if h_ref_sq <= 0
            return;
        end

        ref_direction = P_half \ h_ref;
        alpha = sqrt(max(0, (W_eps - sigma)) / h_ref_sq);
        x_ref_plus = x - e_tran - alpha * ref_direction;

        y_plus(3:4) = x_ref_plus;

        e_plus = x - x_ref_plus;
        eps_plus = e_plus - e_tran;
        Ieps_jump_plus = eps_plus' * p.P * eps_plus;
        y_plus(11) = max(0, Ieps + (Ieps_jump_plus - W_eps));

        p.Ieps_anchor = y_plus(11);
        p.S_eps_ref_anchor = S_eps;
        p.t_last_ref_jump = t;

        p.ref_jump_times(end+1,1) = t;
        p.Ieps_anchor_hist(end+1,1) = y_plus(11);

    case 3
        if W_et <= 0
            return;
        end

        if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
            sigma = min(p.DeltaW_step_tran, W_et);
            beta_fixed = sqrt(max(0, (W_et - sigma)) / W_et);
            e_tran_plus = -beta_fixed * e_tran;
            y_plus(9:10) = e_tran_plus;

            if isfield(p,'fixed_integral_jump_correction') && p.fixed_integral_jump_correction
                Itran_jump_plus = e_tran_plus' * p.P_tran * e_tran_plus;
                y_plus(12) = max(0, Itran + (Itran_jump_plus - W_et));
            end

            p.Itran_anchor = y_plus(12);
            p.tran_jump_times(end+1,1) = t;
            p.Itran_anchor_hist(end+1,1) = y_plus(12);
            return;
        end

        p.s_tran = find_s(p.s_tran, W_et, alpha_series);
        sigma = min(series_element(p.s_tran, alpha_series), W_et);

        P_tran_half = sqrtm(p.P_tran);
        h_tran = p.h_tran(t, e_tran, p.P_tran);
        h_tran_sq = h_tran' * h_tran;
        if h_tran_sq <= 0
            return;
        end

        tran_direction = P_tran_half \ h_tran;
        beta = sqrt(max(0, (W_et - sigma)) / h_tran_sq);
        e_tran_plus = -beta * tran_direction;

        y_plus(9:10) = e_tran_plus;

        Itran_jump_plus = e_tran_plus' * p.P_tran * e_tran_plus;
        y_plus(12) = max(0, Itran + (Itran_jump_plus - W_et));

        p.Itran_anchor = y_plus(12);
        p.S_eps_tran_anchor = S_eps;
        p.S_tran_anchor = S_tran;
        p.t_last_tran_jump = t;

        p.tran_jump_times(end+1,1) = t;
        p.Itran_anchor_hist(end+1,1) = y_plus(12);
end

end
