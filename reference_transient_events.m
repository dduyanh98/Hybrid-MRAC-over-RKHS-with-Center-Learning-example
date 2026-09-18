function [value, isterminal, direction] = reference_transient_events(t, y, p)
% REFERENCE_TRANSIENT_EVENTS  Event surfaces for hybrid jumps.
%   value(1): epoch jump, active for modes 2 and 3
%   value(2): reference jump, Eq. 46, active for mode 3
%   value(3): transient jump, Eq. 48, active for mode 3
%   value(4): kernel-center update, active for mode 3
%   value(5): active kernel-region exit, active for modes 2 and 3

% Epoch event: t = t_last_epoch + T_epoch.
val_epoch = (p.t_last_epoch + p.T_epoch) - t;

% Reference/transient jumps are disabled outside architecture 3.
val_ref  = 1;
val_tran = 1;
val_kernel = 1;
val_box_exit = 1;

% Stop the flow as soon as the state has moved just outside the region
% represented by the current kernel centers.  The small outward tolerance
% avoids repeatedly selecting the old box when the state lies exactly on a
% shared boundary.
if p.architecture_mode >= 2 && isfield(p,'kernel_active_region') && ...
        ~isempty(p.kernel_active_region)
    B = p.kernel_active_region;
    if ~isfield(p,'kernel_box_exit_tol') || isempty(p.kernel_box_exit_tol)
        box_exit_tol = 1e-9;
    else
        box_exit_tol = p.kernel_box_exit_tol;
    end
    theta = y(1);
    theta_dot = y(2);
    margins = [theta - B.theta_L; B.theta_R - theta; ...
               theta_dot - B.d_D; B.d_U - theta_dot];
    val_box_exit = min(margins) + box_exit_tol;
end

if p.architecture_mode == 3 && isfield(p,'epoch_index') && p.epoch_index > 0
    Ieps  = y(11);   % integral of W(epsilon)
    Itran = y(12);   % integral of W(e_tran)
    S_eps = y(p.idx_S_eps);    % sum of epsilon'P epsilon jump variations
    S_tran = y(p.idx_S_tran);  % sum of e_tran'P_tran e_tran variations

    if ~isfield(p,'event_time_eps') || isempty(p.event_time_eps)
        event_time_eps = 1e-9;
    else
        event_time_eps = p.event_time_eps;
    end

    if ~isfield(p,'event_tol') || isempty(p.event_tol)
        event_tol = 1e-10;
    else
        event_tol = p.event_tol;
    end

    if ~isfield(p,'budget_tol') || isempty(p.budget_tol)
        budget_tol = 1e-4;
    else
        budget_tol = p.budget_tol;
    end

    if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
        val_ref  = (p.ref_budget_cum + p.DeltaV_step_ref)  - Ieps;
        val_tran = (p.tran_budget_cum + p.DeltaW_step_tran) - Itran;
        if isfield(p,'idx_kernel_l2') && isfield(p,'kernel_l2_high') && p.kernel_l2_high > 0
            val_kernel = p.kernel_l2_high^2 - y(p.idx_kernel_l2);
        end
        value      = [val_epoch; val_ref; val_tran; val_kernel; val_box_exit];
        isterminal = [1; 1; 1; 1; 1];
        direction  = [-1; -1; -1; -1; -1];
        return;
    end

    t_last_ref = p.t_last_ref_jump;
    t_last_tran = p.t_last_tran_jump;
    use_positive_budget = ~isfield(p,'positive_jump_budget_only') || ...
                          p.positive_jump_budget_only;

    % Eq. (46): trigger when the integral catches the summed Lyapunov budget.
    if use_positive_budget
        ref_budget = max(0, S_eps);
    else
        ref_budget = S_eps;
    end
    ref_event_enabled = ref_budget > budget_tol;

    if ref_event_enabled
        ref_integral = max(0, Ieps);
        remaining_ref = ref_budget - ref_integral;

        ref_time_guard = max(p.t_last_epoch, t_last_ref) ...
                        + event_time_eps - t;

        val_ref = max(ref_time_guard, remaining_ref - event_tol);
    end

    % Eq. (48): trigger when the integral catches the larger Lyapunov budget.
    if use_positive_budget
        tran_budget = max(max(0, S_eps), max(0, S_tran));
    else
        tran_budget = max(S_eps, S_tran);
    end
    tran_event_enabled = tran_budget > budget_tol;

    if tran_event_enabled
        tran_integral = max(0, Itran);
        remaining_tran = tran_budget - tran_integral;

        tran_time_guard = max([p.t_last_epoch,...
                           t_last_ref,...
                           t_last_tran]) ...
                        + event_time_eps - t;

        val_tran = max(tran_time_guard, remaining_tran - event_tol);
    end
end

if p.architecture_mode == 3 && isfield(p,'idx_kernel_l2') && ...
        isfield(p,'kernel_l2_high') && p.kernel_l2_high > 0
    val_kernel = p.kernel_l2_high^2 - y(p.idx_kernel_l2);
end

value      = [val_epoch; val_ref; val_tran; val_kernel; val_box_exit];
isterminal = [1; 1; 1; 1; 1];
direction  = [-1; -1; -1; -1; -1];

end
