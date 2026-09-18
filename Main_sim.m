% ========================================================================
% HYBRID TWO-LAYER MRAC × LEARNING PENDULUM
%   Mode 1: adaptive only          (no kernel, no jumps)
%   Mode 2: trajectory/epoch kernel (no ref/transient jumps, static boxes)
%   Mode 3: hybrid + quadtree       (all jumps + refinement + pruning)
% ========================================================================

addpath(genpath(pwd));
if ~strcmp(getenv('PENDULUM_KEEP_WORKSPACE'), '1')
    clearvars;
end
close all; clc;
rng(0);

if strcmp(getenv('PENDULUM_NO_PLOTS'), '1')
    set(0, 'DefaultFigureVisible', 'off');
else
    set(0, 'DefaultFigureVisible', 'on');
end
set(0, 'DefaultTextInterpreter', 'latex', ...
       'DefaultLegendInterpreter', 'latex', ...
       'DefaultAxesTickLabelInterpreter', 'latex');

% --- Problem setup ------------------------------------------------------
p = hybrid_params();

p.architecture_mode = 3;     % 1, 2, or 3

%% --- Simulation settings -----------------------------------------------
p.Tmax    = 10;
p.T_epoch = 0.1;

% --- Kernel center selection --------------------------------------------
% Active only when p.architecture_mode == 3.
%   1: outward local grid with 2x2/4x4/6x6 centers by L2 error
%   2: original-box grid with 2x2/3x3/5x5 centers by L2 error
%   3: four-center box refinement/expansion by L2 error
p.kernel_mode = 3;
p.kernel_window_epochs = 10;
p.first_two_plots_only = true;

override_json = getenv('PENDULUM_OVERRIDES');
if ~isempty(override_json)
    overrides = jsondecode(override_json);
    override_fields = fieldnames(overrides);
    for k_override = 1:numel(override_fields)
        field = override_fields{k_override};
        p.(field) = overrides.(field);
    end
end

p = apply_kernel_mode_params(p);

if ~isempty(override_json)
    overrides = jsondecode(override_json);
    override_fields = fieldnames(overrides);
    for k_override = 1:numel(override_fields)
        field = override_fields{k_override};
        p.(field) = overrides.(field);
    end
end

lambda_min_Q = min(eig(p.Q));
p.W_fun = @(z) lambda_min_Q * (z(:)' * z(:));
p.P = lyap(p.A_tran', p.Q);
if isfield(p,'Q_tran') && ~isempty(p.Q_tran)
    p.P_tran = lyap(p.A_tran', p.Q_tran);
end

fprintf('\n========== HYBRID LEARNING SIMULATION ==========\n');

% --- Reference signal ---------------------------------------------------
freq = 1.5;
p.desired_output     = @(t) deg2rad(45 + 5*sin(2*pi*freq*t));
p.desired_output_dot = @(t) deg2rad(2*pi*freq*5*cos(2*pi*freq*t));

% --- Initial states -----------------------------------------------------
theta0    = deg2rad(10);
thetadot0 = 0;
x0     = [theta0; thetadot0];
xref0  = [0; 0];
Kx0    = [0; 0];
Kr0    = 0;
Kg0    = [0; 0];
Theta_rkhs0 = zeros(p.kernel_max_centers,1);
Iref0  = 0;
e_tran0 = x0 - xref0;
Ieps0   = 0;
Itran0  = 0;
Theta_hat0 = 0;
S_eps0 = 0;
S_tran0 = 0;
KernelL20 = 0;
X0 = [x0; xref0; Kx0; Kr0; Iref0; e_tran0; Ieps0; Itran0; Theta_hat0; Kg0; Theta_rkhs0; S_eps0; S_tran0; KernelL20];

% --- Kernel containers --------------------------------------------------
Xi    = []; 
alpha = []; 
y_all = [];

% --- Lyapunov budgets / anchors ----------------------------------------
p.ref_budget_cum   = 0;
p.tran_budget_cum  = 0;
p.Ieps_anchor      = 0;
p.Itran_anchor     = 0;
p.t_last_epoch     = 0;
p.t_last_ref_jump  = 0;
p.t_last_tran_jump = 0;

% --- Logging ------------------------------------------------------------
LOG.T = []; LOG.X = [];
Xi_hist      = {};
alpha_hist   = {};
center_times = [];
updates      = 0;
kernel_l2_window = 0;
kernel_window_start = 0;
if isfield(p,'kernel_mode1_grid_sequence') && ~isempty(p.kernel_mode1_grid_sequence)
    mode1_grid_sequence = p.kernel_mode1_grid_sequence(:)';
    mode1_grid_sequence = unique(max(2, round(mode1_grid_sequence)), 'stable');
    mode1_grid_sequence = mode1_grid_sequence(mode1_grid_sequence.^2 <= p.kernel_max_centers);
    if isempty(mode1_grid_sequence)
        mode1_grid_sequence = [2 3 5];
    end
else
    mode1_grid_sequence = [2 3 5];
end
mode1_grid_n = mode1_grid_sequence(1);
mode2_grid_n = 2;
prevMode1RootID = -1;
mode1ActiveRegion = [];
mode1ActiveRootCenter = [];
prevMode2RootID = -1;
mode2ActiveRegion = [];
mode2ActiveRootCenter = [];
mode3ActiveRegion = [];
mode3ActiveRootCenter = [];
prevMode3RootID = -1;
mode3ActiveXi = [];
mode3ActiveY = [];
mode3ActiveKind = '';
p.kernel_active_region = [];

% --- Quadtree initialization (only modes >= 2) --------------------------
Boxes       = [];
Xi_lib      = [];
y_lib       = [];

if p.architecture_mode >= 2
    theta_grid_deg = 0:10:60;         
    thd_grid_deg   = -40:20:160;      
    
    % Convert to rad
    theta_grid_rad = deg2rad(theta_grid_deg);
    thd_grid_rad   = deg2rad(thd_grid_deg);

    % Build library
    [TH, THD] = ndgrid(theta_grid_rad, thd_grid_rad);
    Xi_lib            = [TH(:), THD(:)];
    n_lib             = size(Xi_lib,1);

    y_lib = zeros(n_lib,1);
    for i = 1:n_lib
        y_lib(i) = regressor_vector(Xi_lib(i,:)'); % 
    end

    n_boxes = (numel(theta_grid_rad)-1) * (numel(thd_grid_rad)-1);
    box_template = struct('theta_L', [], 'theta_R', [], ...
                          'd_D', [], 'd_U', [], ...
                          'center', [], 'corner_idx', [], ...
                          'depth', [], 'active', [], 'root_id', []);
    box_list = repmat(box_template, n_boxes, 1);
    box_count = 0;
    for i = 1:(numel(theta_grid_rad)-1)
        for j = 1:(numel(thd_grid_rad)-1)
            box_count = box_count + 1;
            b.theta_L = theta_grid_rad(i);
            b.theta_R = theta_grid_rad(i+1);
            b.d_D     = thd_grid_rad(j);
            b.d_U     = thd_grid_rad(j+1);
            b.center  = [(b.theta_L + b.theta_R)/2, ...
                         (b.d_D     + b.d_U    )/2];

            BL = sub2ind([7 11], i,   j);
            BR = sub2ind([7 11], i+1, j);
            TL = sub2ind([7 11], i,   j+1);
            TR = sub2ind([7 11], i+1, j+1);
            b.corner_idx = [BL; BR; TL; TR];

            b.depth      = 0;
            b.active     = true;
            b.root_id    = box_count;
            box_list(box_count) = b;
        end
    end
    Boxes = box_list;

    % Initialize the kernel dictionary from the box containing the initial
    % trajectory state.  This lets the box-exit event drive center handoffs
    % from t = 0 instead of waiting for the first epoch or hybrid reset.
    x_rep0 = x0.';
    [Xi_box0, y_box0, BoxID0] = select_centers_from_library( ...
        x_rep0, Xi_lib, y_lib, Boxes, p);
    RootID0 = box_root_id(Boxes, BoxID0);

    if p.architecture_mode == 3
        switch p.kernel_mode
            case 1
                rootBox = Boxes(RootID0);
                [Xi, y_all] = quadtree_grid_centers(Boxes(BoxID0), mode2_grid_n);
                mode2ActiveRegion = box_from_centers(Xi, RootID0);
                mode2ActiveRootCenter = rootBox.center;
                prevMode2RootID = RootID0;
                p.kernel_active_region = Boxes(BoxID0);

            case 2
                rootBox = Boxes(RootID0);
                [Xi, y_all] = full_box_grid_centers(rootBox, mode1_grid_n);
                mode1ActiveRegion = box_from_centers(Xi, RootID0);
                mode1ActiveRootCenter = rootBox.center;
                prevMode1RootID = RootID0;
                p.kernel_active_region = rootBox;

            case 3
                Xi = Xi_box0;
                y_all = y_box0;
                mode3ActiveRegion = Boxes(BoxID0);
                mode3ActiveXi = Xi;
                mode3ActiveY = y_all;
                mode3ActiveRootCenter = Boxes(RootID0).center;
                prevMode3RootID = RootID0;
                mode3ActiveKind = 'base';
                p.kernel_active_region = Boxes(BoxID0);

            otherwise
                error('Unknown p.kernel_mode = %d. Use 1, 2, or 3.', p.kernel_mode);
        end
    else
        Xi = Xi_box0;
        y_all = y_box0;
        p.kernel_active_region = Boxes(BoxID0);
    end

    lambda = 1e-9;
    [~,K_gram,~] = kernel_regression(Xi, Xi, y_all, p.sigma, lambda);
    alpha = K_gram \ eye(size(K_gram));

    updates = updates + 1;
    Xi_hist{updates} = Xi;
    alpha_hist{updates} = alpha;
    center_times(end+1) = 0;
end

%% Dynamical systems 
% =======================================================================
% ARCHITECTURE 1 (no jumps, no kernel)
% =======================================================================
if p.architecture_mode == 1
    fprintf('Architecture 1: adaptive only\n');

    opts = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [Tfull, Xfull] = ode45(@(t,x) hybrid_pendulum(t,x,p,[],[]), ...
                           [0 p.Tmax], X0, opts);

    LOG.T = Tfull;
    LOG.X = Xfull;

    fprintf('\n================================================\n');
    fprintf('Simulation done (mode 1).\n');
    fprintf('Kernel updates = 0\n');
    fprintf('================================================\n');

else
% =======================================================================
% ARCHITECTURES 2 & 3 — EVENT- AND TRAJECTORY-BASED KERNEL HANDOFFS
% =======================================================================
    if p.architecture_mode == 2
        fprintf('Architecture 2: trajectory/epoch kernel\n');
    else
        fprintf('Architecture 3: full hybrid + quadtree kernel\n');
    end

    t0 = 0;
    while t0 < p.Tmax

        opts_ev = odeset('RelTol',1e-3, ...
                         'AbsTol',1e-6*ones(length(X0),1), ...
                         'Events', @(t,y) reference_transient_events(t,y,p));

        [Tseg,Xseg,te,~,ie] = ode45(@(t,x) hybrid_pendulum(t,x,p,Xi,alpha), ...
                                    [t0 p.Tmax], X0, opts_ev);

        if p.architecture_mode == 3 && ...
                (~isfield(p,'accumulate_flow_budget') || p.accumulate_flow_budget)
            [Xseg,p] = accumulate_positive_budget_segment(Xseg,p);
        end

        if p.architecture_mode == 3 && isfield(p,'idx_kernel_l2')
            kernel_l2_window = Xseg(end,p.idx_kernel_l2);
        end

        LOG.T = [LOG.T; Tseg];
        LOG.X = [LOG.X; Xseg];

        if isempty(ie) || Tseg(end) >= p.Tmax
            break;
        end

        which   = min(ie);
        t_event = te(end);

        % ---------------------------------------------------------
        % Reset series logic
        % ---------------------------------------------------------
        if p.architecture_mode == 3 && isfield(p,'reset_series_time')
            if p.s_ref > 1e9 || (p.s_ref > 0 ...
                    && t_event - p.t_last_ref_jump > p.reset_series_time)
                p.s_ref = 0;
                fprintf('    Resetting reference series\n');
            end

            if p.s_tran > 1e9 || (p.s_tran > 0 ...
                    && t_event - p.t_last_tran_jump > p.reset_series_time)
                p.s_tran = 0;
                fprintf('    Resetting transient series\n');
            end
        end

        % ---------------------------------------------------------
        % Mode 2 ignores reference/transient and L2 refinement jumps, but
        % still responds immediately when the state leaves its kernel box.
        % ---------------------------------------------------------
        if p.architecture_mode == 2 && any(which == [2 3 4])
            t0 = t_event;
            X0 = Xseg(end,:)';
            continue;
        end

        switch which
            case 1
                fprintf('t = %.2f → Epoch jump\n', t_event);
            case 2
                fprintf('t = %.2f → Reference jump\n', t_event);
            case 3
                fprintf('t = %.2f → Transient jump\n', t_event);
            case 4
                fprintf('t = %.2f -> Kernel-center event\n', t_event);
            case 5
                fprintf('t = %.2f -> Kernel-box exit event\n', t_event);
        end

        % ---------------------------------------------------------
        % SAVE PRE-JUMP STATE
        % ---------------------------------------------------------
        y_minus = Xseg(end,:)';
        p_minus = p;

        % ---------------------------------------------------------
        % APPLY JUMP
        % ---------------------------------------------------------
        kernel_event_only = any(which == [4 5]);
        if kernel_event_only
            X_plus = y_minus;
        elseif which == 1 || p.architecture_mode == 3
            [X_plus,p] = apply_jump(which, t_event, y_minus, p);
        else
            X_plus = y_minus;
        end

        % ---------------------------------------------------------
        % JUMP-BASED LYAPUNOV BUDGET UPDATE
        % ---------------------------------------------------------
        if ~kernel_event_only
        dVeps = lyapunov_eps_value(X_plus,p) ...
            - lyapunov_eps_value(y_minus,p_minus);

        dVtran = lyapunov_tran_value(X_plus,p) ...
            - lyapunov_tran_value(y_minus,p_minus);

        if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
            if ~isfield(p,'fixed_budget_law') || isempty(p.fixed_budget_law)
                fixed_budget_law = 'constant';
            else
                fixed_budget_law = p.fixed_budget_law;
            end

            switch fixed_budget_law
                case 'constant'
                    dVeps_budget = p.DeltaV_step_ref;
                    dVtran_budget = p.DeltaW_step_tran;
                case {'floor', 'floor_actual'}
                    dVeps_budget = max(p.DeltaV_step_ref, dVeps);
                    dVtran_budget = max(p.DeltaW_step_tran, dVtran);
                case {'positive', 'positive_actual'}
                    dVeps_budget = max(0, dVeps);
                    dVtran_budget = max(0, dVtran);
                otherwise
                    error('Unknown p.fixed_budget_law = %s', fixed_budget_law);
            end

            if which == 2
                p.ref_budget_cum = p.ref_budget_cum + dVeps_budget;
                p.DeltaV_hist(end+1,1) = dVeps_budget;
            elseif which == 3
                p.tran_budget_cum = p.tran_budget_cum + dVtran_budget;
                p.DeltaW_hist(end+1,1) = dVtran_budget;
            end
        else
            if ~isfield(p,'positive_jump_budget_only') || p.positive_jump_budget_only
                if ~isfield(p,'min_jump_budget_eps') || isempty(p.min_jump_budget_eps)
                    min_jump_budget_eps = 0;
                else
                    min_jump_budget_eps = p.min_jump_budget_eps;
                end

                if ~isfield(p,'min_jump_budget_tran') || isempty(p.min_jump_budget_tran)
                    min_jump_budget_tran = 0;
                else
                    min_jump_budget_tran = p.min_jump_budget_tran;
                end

                dVeps_budget = max(min_jump_budget_eps, dVeps);
                dVtran_budget = max(min_jump_budget_tran, dVtran);
            else
                dVeps_budget = dVeps;
                dVtran_budget = dVtran;
            end

            if dVeps_budget ~= 0
                X_plus(p.idx_S_eps) = X_plus(p.idx_S_eps) + dVeps_budget;
            end

            if dVtran_budget ~= 0
                X_plus(p.idx_S_tran) = X_plus(p.idx_S_tran) + dVtran_budget;
            end

            p.DeltaV_hist(end+1,1) = dVeps_budget;
            p.DeltaW_hist(end+1,1) = dVtran_budget;
            p.budget_jump_times(end+1,1) = t_event;
        end
        end
        % ----------------- KERNEL UPDATE (modes 2 and 3) --------------------
        % Now: update kernel after *any* jump (epoch, ref, transient) in mode >= 2
        if p.architecture_mode >= 2


            x_rep = Xseg(end,1:2);
            if which ~= 5 && isfield(p,'kernel_rep_point') && ...
                    strcmp(p.kernel_rep_point,'max_tracking_error')
                local_error_signal = sum((Xseg(:,1:2) - Xseg(:,3:4)).^2, 2);
                [~, idx_bad] = max(local_error_signal);
                x_rep = Xseg(idx_bad,1:2);
            end
            kernel_l2 = sqrt(max(kernel_l2_window, 0));

            [~, ~, BoxID0] = select_centers_from_library(x_rep, Xi_lib, y_lib, Boxes, p);
            RootID0 = box_root_id(Boxes, BoxID0);

            if p.architecture_mode == 3
                switch p.kernel_mode
                    case 2
                        rootChanged = RootID0 ~= prevMode1RootID;
                        rootBox = Boxes(RootID0);
                        if rootChanged && ~isempty(mode1ActiveRegion)
                            [Xi_new, y_new] = translate_kernel_centers( ...
                                Xi, mode1ActiveRootCenter, rootBox.center);
                            mode1ActiveRegion = box_from_centers(Xi_new, RootID0);
                        else
                            if ~rootChanged && kernel_l2 > p.kernel_l2_high
                                idx_grid = find(mode1_grid_sequence > mode1_grid_n, 1, 'first');
                                if ~isempty(idx_grid)
                                    mode1_grid_n = mode1_grid_sequence(idx_grid);
                                end
                            elseif ~rootChanged && kernel_l2 < p.kernel_l2_low
                                idx_grid = find(mode1_grid_sequence < mode1_grid_n, 1, 'last');
                                if ~isempty(idx_grid)
                                    mode1_grid_n = mode1_grid_sequence(idx_grid);
                                end
                            end
                            [Xi_new, y_new] = full_box_grid_centers(rootBox, mode1_grid_n);
                            mode1ActiveRegion = box_from_centers(Xi_new, RootID0);
                        end
                        prevMode1RootID = RootID0;
                        mode1ActiveRootCenter = rootBox.center;
                        BoxID = RootID0;
                        fprintf('    Kernel mode 2: root box %d, grid=%dx%d, centers=%d, L2=%.4g\n', ...
                                RootID0, mode1_grid_n, mode1_grid_n, size(Xi_new,1), kernel_l2);

                    case 1
                        rootChanged = RootID0 ~= prevMode2RootID;
                        rootBox = Boxes(RootID0);

                        if rootChanged && ~isempty(mode2ActiveRegion)
                            [Xi_new, y_new] = translate_kernel_centers( ...
                                Xi, mode2ActiveRootCenter, rootBox.center);
                        else
                            if ~rootChanged && kernel_l2 > p.kernel_l2_high
                                mode2_grid_n = min(mode2_grid_n + 2, p.kernel_mode2_max_grid_n);
                            elseif ~rootChanged && kernel_l2 < p.kernel_l2_low
                                mode2_grid_n = max(mode2_grid_n - 2, p.kernel_mode2_min_grid_n);
                            end
                            [Xi_new, y_new] = quadtree_grid_centers(Boxes(BoxID0), mode2_grid_n);
                        end

                        mode2ActiveRegion = box_from_centers(Xi_new, RootID0);
                        mode2ActiveRootCenter = rootBox.center;
                        prevMode2RootID = RootID0;
                        BoxID = BoxID0;
                        fprintf('    Kernel mode 1: grid=%dx%d, centers=%d, L2=%.4g\n', ...
                                mode2_grid_n, mode2_grid_n, size(Xi_new,1), kernel_l2);

                    case 3
                        [Xi_box, y_box, BoxID] = select_centers_from_library( ...
                            x_rep, Xi_lib, y_lib, Boxes, p);
                        BoxID_leaf = BoxID;
                        PruneBoxID = BoxID;
                        RootID = box_root_id(Boxes, BoxID);
                        rootChanged = RootID ~= prevMode3RootID;

                        if rootChanged && ~isempty(mode3ActiveXi)
                            rootBox = Boxes(RootID);
                            [Xi_new, y_new] = translate_kernel_centers( ...
                                mode3ActiveXi, mode3ActiveRootCenter, rootBox.center);
                            BoxID = RootID;
                            PruneBoxID = BoxID_leaf;
                            mode3ActiveRegion = box_from_centers(Xi_new, RootID);
                            mode3ActiveXi = Xi_new;
                            mode3ActiveY = y_new;
                            mode3ActiveRootCenter = rootBox.center;
                            prevMode3RootID = RootID;
                            fprintf('    Kernel mode 3: shifted active centers to root box %d, L2=%.4g\n', ...
                                    RootID, kernel_l2);
                        elseif rootChanged
                            rootBox = Boxes(RootID);
                            Xi_new = Xi_lib(rootBox.corner_idx, :);
                            y_new = y_lib(rootBox.corner_idx);
                            BoxID = RootID;
                            PruneBoxID = BoxID_leaf;
                            mode3ActiveRegion = rootBox;
                            mode3ActiveXi = Xi_new;
                            mode3ActiveY = y_new;
                            mode3ActiveRootCenter = rootBox.center;
                            prevMode3RootID = RootID;
                            mode3ActiveKind = 'base';
                            fprintf('    Kernel mode 3: new root box %d, using original centers, L2=%.4g\n', ...
                                    RootID, kernel_l2);
                        elseif strcmp(mode3ActiveKind, 'expanded')
                            Xi_new = mode3ActiveXi;
                            y_new = mode3ActiveY;
                            fprintf('    Kernel mode 3: keeping expanded box centers, L2=%.4g\n', kernel_l2);
                        elseif kernel_l2 > p.kernel_l2_high && Boxes(BoxID).depth < p.max_depth
                            fprintf('    Kernel mode 3: refining box %d (depth=%d), L2=%.4g\n', ...
                                    BoxID, Boxes(BoxID).depth, kernel_l2);
                            [Boxes, Xi_lib, y_lib] = refine_box(Boxes, BoxID, Xi_lib, y_lib, p);
                            [Xi_new, y_new, BoxID] = select_centers_from_library( ...
                                x_rep, Xi_lib, y_lib, Boxes, p);
                            mode3ActiveRegion = Boxes(BoxID);
                            mode3ActiveXi = Xi_new;
                            mode3ActiveY = y_new;
                            mode3ActiveRootCenter = Boxes(RootID).center;
                            prevMode3RootID = RootID;
                            mode3ActiveKind = 'refined';
                        elseif kernel_l2 < p.kernel_l2_low
                            [Xi_new, y_new] = expand_box_centers(x_rep, Xi_box, Xi_lib, y_lib);
                            mode3ActiveRegion = box_from_centers(Xi_new, RootID);
                            mode3ActiveXi = Xi_new;
                            mode3ActiveY = y_new;
                            mode3ActiveRootCenter = Boxes(RootID).center;
                            prevMode3RootID = RootID;
                            mode3ActiveKind = 'expanded';
                            fprintf('    Kernel mode 3: expanded box centers, L2=%.4g\n', kernel_l2);
                        else
                            Xi_new = Xi_box;
                            y_new = y_box;
                            mode3ActiveRegion = Boxes(BoxID);
                            mode3ActiveXi = Xi_new;
                            mode3ActiveY = y_new;
                            mode3ActiveRootCenter = Boxes(RootID).center;
                            prevMode3RootID = RootID;
                            mode3ActiveKind = 'base';
                            fprintf('    Kernel mode 3: box %d (depth=%d), L2=%.4g\n', ...
                                    BoxID, Boxes(BoxID).depth, kernel_l2);
                        end

                        Boxes = prune_boxes(x_rep, Boxes, t_event, p, PruneBoxID);

                    otherwise
                        error('Unknown p.kernel_mode = %d. Use 1, 2, or 3.', p.kernel_mode);
                end

            else
                [Xi_new, y_new, BoxID] = select_centers_from_library( ...
                    x_rep, Xi_lib, y_lib, Boxes, p);
                fprintf('    Using box %d (depth=%d)\n', ...
                        BoxID, Boxes(BoxID).depth);
            end

            % Publish the current spatial box so the next ODE segment
            % terminates immediately upon departure.  Modes 1 and 2 can use
            % center stencils wider than the box; those wider stencils must
            % not delay the handoff to the neighboring box.
            if p.architecture_mode == 2
                p.kernel_active_region = Boxes(BoxID);
            else
                switch p.kernel_mode
                    case 1
                        p.kernel_active_region = Boxes(BoxID0);
                    case 2
                        p.kernel_active_region = rootBox;
                    case 3
                        p.kernel_active_region = Boxes(BoxID);
                end
            end

            Xi    = Xi_new;
            y_all = y_new;
            lambda = 1e-9;
            [~,K_gram,~] = kernel_regression(Xi, Xi, y_all, p.sigma, lambda);
            alpha = K_gram \ eye(size(K_gram));

            updates             = updates + 1;
            Xi_hist{updates}    = Xi;
            alpha_hist{updates} = alpha;
            center_times(end+1) = t_event;
        end
        % ===================================================================

        if p.architecture_mode == 3 && isfield(p,'idx_kernel_l2')
            X_plus(p.idx_kernel_l2) = 0;
            kernel_l2_window = 0;
            kernel_window_start = t_event;
        elseif p.architecture_mode == 3 && ...
                t_event - kernel_window_start >= p.kernel_window_epochs*p.T_epoch - 1e-12
            kernel_window_start = t_event;
            kernel_l2_window = 0;
        end

        if which == 1
            p.t_last_epoch = t_event;
        end

        t0 = t_event;
        X0 = X_plus;
        if t0 >= p.Tmax, break; end
    end

    fprintf('\n================================================\n');
    fprintf('Simulation done (mode %d).\n', p.architecture_mode);
    fprintf('Total kernel updates = %d\n', updates);
    fprintf('================================================\n');

end

% =======================================================================
% POST–PROCESSING & PLOTS  (CLEAN + CONSISTENT)
% =======================================================================
T = LOG.T; 
X = LOG.X;
N = length(T);

% State extraction
x      = X(:,1:2);
xref   = X(:,3:4);
Kx     = X(:,5:6);
Kr     = X(:,7);
Iref   = X(:,8);
e_tran = X(:,9:10);
Ieps   = X(:,11);       % ∫ W(ε)
Itran  = X(:,12);       % ∫ W(e_tran)
Ex_f_hat_state = X(:,13);
Kg     = X(:,14:15);
Theta_rkhs = X(:,p.idx_theta_rkhs);
S_eps  = X(:,p.idx_S_eps);
S_tran = X(:,p.idx_S_tran);

% Errors
e   = x - xref;
eps = e - e_tran;

theta     = x(:,1);
theta_dot = x(:,2);
theta_ref     = xref(:,1);
theta_ref_dot = xref(:,2);

% Plot trajectories using the last sample at duplicate hybrid-event times.
% This shows post-jump values and avoids visual spikes from pre/post jump pairs.
[T_plot, idx_plot_last] = unique(T, 'last');
[T_plot, idx_plot_order] = sort(T_plot);
idx_plot_last = idx_plot_last(idx_plot_order);

theta_plot         = theta(idx_plot_last);
theta_dot_plot     = theta_dot(idx_plot_last);
theta_ref_plot     = theta_ref(idx_plot_last);
theta_ref_dot_plot = theta_ref_dot(idx_plot_last);
theta_des_plot     = arrayfun(p.desired_output, T_plot);
theta_des_dot_plot = arrayfun(p.desired_output_dot, T_plot);

if strcmp(getenv('PENDULUM_SKIP_PLOT_BUILD'), '1')
    return;
end

mode1_center_change_times = [];
mode1_center_change_counts = [];
mode1_center_count_values = [];
mode1_center_count_colors = [];
show_center_change_lines = isfield(p, 'show_center_change_lines') && p.show_center_change_lines;
if show_center_change_lines && p.architecture_mode == 3 && p.kernel_mode == 2 && numel(Xi_hist) >= 2
    mode1_center_counts = cellfun(@(Xi_k) size(Xi_k, 1), Xi_hist(:));
    mode1_change_idx = find(diff(mode1_center_counts) ~= 0) + 1;
    if ~isempty(mode1_change_idx)
        mode1_center_change_times = center_times(mode1_change_idx);
        mode1_center_change_counts = mode1_center_counts(mode1_change_idx);
        mode1_center_count_values = unique(mode1_center_counts, 'stable');
        mode1_center_count_colors = lines(numel(mode1_center_count_values));
    end
end

% =======================================================================
% BASIC TRAJECTORY PLOTS
% =======================================================================
switch p.kernel_mode
    case 1
        center_mode_label = 'kernel center mode 1';
    case 2
        center_mode_label = 'kernel center mode 2';
    case 3
        center_mode_label = 'kernel center mode 3';
    otherwise
        center_mode_label = sprintf('kernel center mode %d', p.kernel_mode);
end
center_title_suffix = sprintf('Center: %s', center_mode_label);
legend_font_size = 15;
legend_box_scale = 1.12;

figure('Color','w');
subplot(2,1,1);
h_theta = plot(T_plot, rad2deg(theta_plot),'k','LineWidth',1.5); hold on;
h_theta_ref = plot(T_plot, rad2deg(theta_ref_plot),'r--','LineWidth',1.5);
% h_theta_des = plot(T_plot, rad2deg(theta_des_plot),'b:','LineWidth',1.8);
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\theta$ [deg]', 'Interpreter', 'latex');
h_legend = legend([h_theta h_theta_ref %h_theta_des
    ], ...
    {'$\theta$', '$\theta_{\mathrm{ref}}$'},'Location','best','Interpreter','latex','FontSize',legend_font_size);
enlarge_legend_box(h_legend, legend_box_scale);
grid on;
title(sprintf('Tracking of $\\theta$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

subplot(2,1,2);
h_theta_dot = plot(T_plot, rad2deg(theta_dot_plot),'k','LineWidth',1.5); hold on;
h_theta_ref_dot = plot(T_plot, rad2deg(theta_ref_dot_plot),'r--','LineWidth',1.5);
% h_theta_des_dot = plot(T_plot, rad2deg(theta_des_dot_plot),'b:','LineWidth',1.8);
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [deg/s]', 'Interpreter', 'latex');
h_legend = legend([h_theta_dot h_theta_ref_dot], ...
    {'$\dot{\theta}$', '$\dot{\theta}_{\mathrm{ref}}$'}, ...
    'Location','best','Interpreter','latex','FontSize',legend_font_size);
enlarge_legend_box(h_legend, legend_box_scale);
grid on;
title(sprintf('Tracking of $\\dot{\\theta}$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

figure('Color','w');
subplot(2,1,1);
h_theta_err = plot(T_plot, rad2deg(theta_plot - theta_ref_plot),'k','LineWidth',1.5); hold on;
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_\theta$ [deg]', 'Interpreter', 'latex');
h_legend = legend(h_theta_err, {'$\theta - \theta_{\mathrm{ref}}$'}, ...
    'Location','best','Interpreter','latex','FontSize',legend_font_size);
enlarge_legend_box(h_legend, legend_box_scale);
grid on;
title(sprintf('Error in $\\theta$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

subplot(2,1,2);
h_theta_dot_err = plot(T_plot, rad2deg(theta_dot_plot - theta_ref_dot_plot),'k','LineWidth',1.5); hold on;
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_{\dot{\theta}}$ [deg/s]', 'Interpreter', 'latex');
h_legend = legend(h_theta_dot_err, {'$\dot{\theta} - \dot{\theta}_{\mathrm{ref}}$'}, ...
    'Location','best','Interpreter','latex','FontSize',legend_font_size);
enlarge_legend_box(h_legend, legend_box_scale);
grid on;
title(sprintf('Error in $\\dot{\\theta}$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

figure('Color','w');
subplot(2,1,1);
h_theta_des_err = plot(T_plot, rad2deg(theta_plot - theta_des_plot),'k','LineWidth',1.5); hold on;
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_{\theta,\mathrm{des}}$ [deg]', 'Interpreter', 'latex');
legend(h_theta_des_err, {'$\theta - \theta_{\mathrm{des}}$'}, ...
    'Location','best','Interpreter','latex','FontSize',legend_font_size);
grid on;
title(sprintf('Desired error in $\\theta$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

subplot(2,1,2);
h_theta_des_dot_err = plot(T_plot, rad2deg(theta_dot_plot - theta_des_dot_plot),'k','LineWidth',1.5); hold on;
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('$e_{\dot{\theta},\mathrm{des}}$ [deg/s]', 'Interpreter', 'latex');
legend(h_theta_des_dot_err, {'$\dot{\theta} - \dot{\theta}_{\mathrm{des}}$'}, ...
    'Location','best','Interpreter','latex','FontSize',legend_font_size);
grid on;
title(sprintf('Desired error in $\\dot{\\theta}$ (%s)', center_title_suffix), ...
    'Interpreter', 'latex');

if isfield(p, 'first_two_plots_only') && p.first_two_plots_only
    return;
end

% =======================================================================
% CONTROL DECOMPOSITION
% =======================================================================
u          = zeros(N,1);
u_baseline = zeros(N,1);
u_adaptive = zeros(N,1);

for i = 1:N
    t = T(i);
    yd   = p.desired_output(t);
    yd_d = p.desired_output_dot(t);

    xi    = x(i,:)';
    Kx_i  = Kx(i,:)';
    Kr_i  = Kr(i);
    Kg_i  = Kg(i,:)';
    Theta_rkhs_i = Theta_rkhs(i,:)';
    Iref_i= Iref(i);
    xr_i = xref(i,:)';
    e_i = xi - xr_i;

    u_bl = p.m_est*p.g*sin(xi(1)) + p.k_est*p.l_est*xi(2) ...
          - p.m_est*p.l_est*(p.kp*(xi(1)-yd) + p.kd*(xi(2)-yd_d));

    % kernel snapshot if available
    if ~isempty(center_times)
        idx = find(center_times <= t, 1, 'last');
    else
        idx = [];
    end

    if ~isempty(idx) && idx <= numel(Xi_hist)
        [~, K_aprx_i] = kernel_prediction(xi', Xi_hist{idx}, zeros(size(alpha_hist{idx},1),1), p.sigma, 1);
        Phi_rkhs_i = K_aprx_i(:);
        n_rkhs_i = min(numel(Theta_rkhs_i), numel(Phi_rkhs_i));
        Ex_f_hat_i = Theta_rkhs_i(1:n_rkhs_i).' * Phi_rkhs_i(1:n_rkhs_i);
    else
        Ex_f_hat_i = Ex_f_hat_state(i);
    end

    r_i = p.kp*yd - p.ki*Iref_i;
    u_ad = Kx_i'*xi + Kr_i*r_i - Ex_f_hat_i;
    if p.architecture_mode == 3
        u_ad = u_ad + Kg_i' * e_i;
    end

    u(i)          = u_bl + u_ad;
    u_baseline(i) = u_bl;
    u_adaptive(i) = u_ad;
end

figure('Color','w');
plot(T,u,'k','LineWidth',1.8); hold on;
plot(T,u_baseline,'r--','LineWidth',1.3);
plot(T,u_adaptive,'b:','LineWidth',1.8);
ylabel('$u$', 'Interpreter', 'latex');
xlabel('$t$ [s]', 'Interpreter', 'latex');
legend({'$u_{\mathrm{total}}$', '$u_{\mathrm{baseline}}$', ...
    '$u_{\mathrm{adaptive}}$'}, ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', legend_font_size);
grid on;
title('Control Effort Decomposition', 'Interpreter', 'latex');

% =======================================================================
% PHASE PORTRAIT + CENTER EVOLUTION
% =======================================================================
figure('Color','w');
plot(rad2deg(theta), rad2deg(theta_dot),'k','LineWidth',1.5); hold on;
plot(rad2deg(theta_ref), rad2deg(theta_ref_dot),'r--','LineWidth',1.2);
xlabel('$\theta$ [deg]', 'Interpreter', 'latex');
ylabel('$\dot{\theta}$ [deg/s]', 'Interpreter', 'latex');
legend({'Actual', 'Reference'}, ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', legend_font_size);
grid on;
title('Phase Portrait', 'Interpreter', 'latex');

% ---- Overlay center locations (if any) --------------------------------
if p.architecture_mode >= 2 && ~isempty(Xi)
    plot(rad2deg(Xi(:,1)), rad2deg(Xi(:,2)), 'bo', ...
         'MarkerFaceColor','c', 'MarkerSize',7, ...
         'DisplayName','Centers');
    for j = 1:size(Xi,1)
        text(rad2deg(Xi(j,1)) + 0.5, rad2deg(Xi(j,2)), sprintf('%d', j), ...
            'FontSize',8, 'Color','b', 'Interpreter', 'latex');
    end
    legend('show','Location','best','Interpreter','latex', ...
        'FontSize', legend_font_size);
end

% ---- Overlay evolution of 4 centers used per update --------------------
if ~isempty(Xi_hist)
    cmap = lines(numel(Xi_hist));
    for k = 1:numel(Xi_hist)
        Xi_k = Xi_hist{k};   % 4×2
        plot(rad2deg(Xi_k(:,1)), rad2deg(Xi_k(:,2)), 'o', ...
             'MarkerSize',6,'Color',cmap(k,:));

        % place label at average of 4 centers
        xc = mean(rad2deg(Xi_k(:,1)));
        yc = mean(rad2deg(Xi_k(:,2)));
        text(xc, yc, sprintf('%d',k), 'FontSize',7, ...
            'Color',cmap(k,:), 'Interpreter', 'latex');
    end
end


% =======================================================================
% OVERLAY CENTER TIMES (if any)
% =======================================================================
% if ~isempty(center_times)
%     figs = [1 2 3];
%     for f = figs
%         figure(f); hold on;
%         for s = 1:length(center_times)
%             xline(center_times(s),'k--','LineWidth',0.8,'HandleVisibility','off');
%         end
%         xline(NaN,'k--','DisplayName','Kernel Update'); % legend handle
%     end
% end

% =======================================================================
% COMPUTE LYAPUNOV QUANTITIES (CLEAN VERSION)
% =======================================================================
Gamma = blkdiag(p.Gamma_x, p.Gamma_r);
if p.architecture_mode == 3
    Gamma = blkdiag(Gamma, p.Gamma_g);
end

theta_vec = [Kx, Kr];  % (N × 3)
if p.architecture_mode == 3
    theta_vec = [theta_vec, Kg];
end

V_eps        = zeros(N,1);  % εᵀ P ε
V_e          = zeros(N,1);  % eᵀ P e
V_theta      = zeros(N,1);  % θ-parameter error term
V_main       = zeros(N,1);  % V = εᵀ P ε + θᵀ Γ^{-1} θ
V_ref_monitor  = zeros(N,1); 
V_tran_monitor = zeros(N,1);

if isfield(p,'ref_jump_times') && ~isempty(p.ref_jump_times)
    ref_times   = p.ref_jump_times(:);
else
    ref_times   = [];
end

if isfield(p,'Ieps_anchor_hist') && ~isempty(p.Ieps_anchor_hist)
    ref_anchors = p.Ieps_anchor_hist(:);
else
    ref_anchors = [];
end

if isfield(p,'tran_jump_times') && ~isempty(p.tran_jump_times)
    tran_times  = p.tran_jump_times(:);
else
    tran_times  = [];
end

if isfield(p,'Itran_anchor_hist') && ~isempty(p.Itran_anchor_hist)
    tran_anchors = p.Itran_anchor_hist(:);
else
    tran_anchors = [];
end


ref_idx = 0; 
tran_idx = 0;
ref_anchor = 0;
tran_anchor = 0;
tol = 1e-9;

for i = 1:N
    eps_i    = eps(i,:)';
    e_tran_i = e_tran(i,:)';
    theta_err= theta_vec(i,:)';

    V_eps(i)   = eps_i' * p.P * eps_i;
    V_e(i)     = e(i,:) * p.P * e(i,:)';
    V_theta(i) = theta_err' * (Gamma \ theta_err);
    V_main(i)  = lyapunov_main_value(X(i,:)', p);

    t_i = T(i);

    % update reference anchor
    if ~isempty(ref_times)
        while ref_idx+1 <= numel(ref_times) && t_i >= ref_times(ref_idx+1)-tol
            ref_idx = ref_idx+1;
            ref_anchor = ref_anchors(ref_idx);
        end
    end

    % update transient anchor
    if ~isempty(tran_times)
        while tran_idx+1 <= numel(tran_times) && t_i >= tran_times(tran_idx+1)-tol
            tran_idx = tran_idx+1;
            tran_anchor = tran_anchors(tran_idx);
        end
    end

    V_ref_monitor(i)  = Ieps(i)  - ref_anchor; 
    V_tran_monitor(i) = Itran(i) - tran_anchor;
end

% =======================================================================
% FIGURE 5 — MAIN LYAPUNOV FUNCTION
% =======================================================================
SIM_METRICS = struct();
SIM_METRICS.tracking_l2 = sqrt(trapz(T, sum(e.^2, 2)));
SIM_METRICS.theta_rmse_deg = sqrt(mean(rad2deg(theta - theta_ref).^2));
SIM_METRICS.theta_dot_rmse_deg = sqrt(mean(rad2deg(theta_dot - theta_ref_dot).^2));
SIM_METRICS.max_abs_theta_err_deg = max(abs(rad2deg(theta - theta_ref)));
SIM_METRICS.control_rms = sqrt(mean(u.^2));
SIM_METRICS.control_max_abs = max(abs(u));
SIM_METRICS.final_V_main = V_main(end);
SIM_METRICS.final_Ve = V_e(end);
SIM_METRICS.final_Veps = V_eps(end);
SIM_METRICS.Ve_l2 = sqrt(trapz(T, V_e.^2));
SIM_METRICS.Veps_l2 = sqrt(trapz(T, V_eps.^2));
SIM_METRICS.kernel_updates = updates;

if strcmp(getenv('PENDULUM_NO_PLOTS'), '1')
    close all;
    return;
end

figure('Color','w');
hold on;
semilogy(T, V_main,'k','LineWidth',1.5);
semilogy(T, V_e,'b','LineWidth',1.4);
semilogy(T, V_eps,'r--','LineWidth',1.4);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('Lyapunov value', 'Interpreter', 'latex');
title('Main Lyapunov Function', 'Interpreter', 'latex');
legend({'$V_{\mathrm{main}}$', '$V_e = e^T P e$', ...
    '$V_{\epsilon} = \epsilon^T P \epsilon$'}, ...
    'Location','best','Interpreter','latex', 'FontSize', legend_font_size);
grid on;

theta_error_l2_cum = sqrt(cumtrapz(T, (theta - theta_ref).^2));

figure('Color','w'); hold on;
plot(T, theta_error_l2_cum,'b','LineWidth',1.5, ...
    'DisplayName', '$\|\theta-\theta_{\mathrm{ref}}\|_{L_2(0,t)}$');
add_mode1_center_change_lines(mode1_center_change_times, ...
    mode1_center_change_counts, mode1_center_count_values, ...
    mode1_center_count_colors);
xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('Cumulative position-error $L_2$ norm [rad$\sqrt{\mathrm{s}}$]', ...
    'Interpreter', 'latex');
title('Position Error $L_2$ Norm', 'Interpreter', 'latex');
legend('show','Location','best','Interpreter','latex', ...
    'FontSize', legend_font_size);
grid on;

% % FIGURE 6 — Eq (46): global integral from t0
% figure('Color','w'); hold on;
% plot(T, Ieps, 'k-', 'LineWidth', 1.5);   % <--- use Ieps directly
% 
% 
% hasDV = isfield(p,'DeltaV_hist') && ~isempty(p.DeltaV_hist);
% if hasDV
%     tJ = p.ref_jump_times(:);
%     DV = p.DeltaV_hist(:);
%     RHS = [0; cumsum(DV)] + p.DeltaV_step_ref;
%     stairs([tJ; T(end)], RHS,'r','LineWidth',1.6);
%     for t0 = tJ'; xline(t0,'r--'); end
% else
%     stairs([T(1);T(end)], [p.DeltaV_step_ref; p.DeltaV_step_ref],'r','LineWidth',1.6);
% end
% 
% xlabel('Time [s]');
% ylabel('Eq. 46 terms');
% legend('LHS = ∫W(ε)', 'RHS = ΣΔV + ΔV_{step}');
% title('Reference Reset Trigger Condition');
% grid on;

% =======================================================================
% FIGURE 6 - Event quantities for reference jumps (Eq. 46)
% =======================================================================

figure('Color','w'); hold on;

if isfield(p,'ref_jump_times') && ~isempty(p.ref_jump_times)
    ref_times = p.ref_jump_times(:);
else
    ref_times = [];
end

lhs_ref_event = Ieps;
if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
    if isfield(p,'DeltaV_hist') && ~isempty(p.DeltaV_hist) && ~isempty(ref_times)
        rhs_ref_event = zeros(N,1);
        ref_idx = 0;
        DV_ref = p.DeltaV_hist(:);

        for i = 1:N
            t_i = T(i);
            while ref_idx < numel(ref_times) && ...
                    ref_idx < numel(DV_ref) && ...
                    t_i >= ref_times(ref_idx+1) - 1e-12
                ref_idx = ref_idx + 1;
            end
            rhs_ref_event(i) = sum(DV_ref(1:ref_idx)) + p.DeltaV_step_ref;
        end
    else
        rhs_ref_event = p.DeltaV_step_ref * ones(N,1);
    end
else
    rhs_ref_event = S_eps;
end

plot(T, lhs_ref_event, 'k-', 'LineWidth', 1.8);
stairs(T, rhs_ref_event, 'r-', 'LineWidth', 1.6);

for t0 = ref_times'
    xline(t0, 'r--', 'LineWidth', 1.1);
end

xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('Eq. 46 Quantities', 'Interpreter', 'latex');
legend({'LHS $= I_{\epsilon}$', 'RHS $=$ active reference budget'}, ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', legend_font_size);
title('Reference Jump', 'Interpreter', 'latex');
grid on;


% =======================================================================
% FIGURE 7 - Event quantities for transient jumps (Eq. 48)
% =======================================================================

figure('Color','w'); hold on;

if isfield(p,'tran_jump_times') && ~isempty(p.tran_jump_times)
    tran_times = p.tran_jump_times(:);
else
    tran_times = [];
end

lhs_tran_event = Itran;
if isfield(p,'hybrid_budget_mode') && strcmp(p.hybrid_budget_mode,'fixed')
    if isfield(p,'DeltaW_hist') && ~isempty(p.DeltaW_hist) && ~isempty(tran_times)
        rhs_tran_event = zeros(N,1);
        tran_idx = 0;
        DW_tran = p.DeltaW_hist(:);

        for i = 1:N
            t_i = T(i);
            while tran_idx < numel(tran_times) && ...
                    tran_idx < numel(DW_tran) && ...
                    t_i >= tran_times(tran_idx+1) - 1e-12
                tran_idx = tran_idx + 1;
            end
            rhs_tran_event(i) = sum(DW_tran(1:tran_idx)) + p.DeltaW_step_tran;
        end
    else
        rhs_tran_event = p.DeltaW_step_tran * ones(N,1);
    end
else
    rhs_tran_event = max(S_eps, S_tran);
end

plot(T, lhs_tran_event, 'k-', 'LineWidth', 1.8);
stairs(T, rhs_tran_event, 'b-', 'LineWidth', 1.6);

for t0 = tran_times'
    xline(t0, 'b--', 'LineWidth', 1.1);
end

xlabel('$t$ [s]', 'Interpreter', 'latex');
ylabel('Eq. 48 Quantities', 'Interpreter', 'latex');
legend({'LHS $= I_{\mathrm{tran}}$', 'RHS $=$ active transient budget'}, ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', legend_font_size);
title('Transient Jump', 'Interpreter', 'latex');
grid on;

function add_mode1_center_change_lines(change_times, change_counts, count_values, count_colors)
if isempty(change_times)
    return;
end

used_count_legend = false(size(count_values));
for kk = 1:numel(change_times)
    count_idx = find(count_values == change_counts(kk), 1);
    if isempty(count_idx)
        continue;
    end

    line_name = sprintf('$%d\\,\\mathrm{centers}$', change_counts(kk));
    if ~used_count_legend(count_idx)
        xline(change_times(kk), ':', line_name, ...
            'Color', count_colors(count_idx,:), ...
            'LineWidth', 1.2, ...
            'Interpreter', 'latex', ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', ...
            'HandleVisibility', 'on', ...
            'DisplayName', line_name);
        used_count_legend(count_idx) = true;
    else
        xline(change_times(kk), ':', ...
            'Color', count_colors(count_idx,:), ...
            'LineWidth', 1.2, ...
            'HandleVisibility', 'off');
    end
end
end

function enlarge_legend_box(lgd, scale)
if isempty(lgd) || ~isvalid(lgd)
    return;
end

drawnow;
old_units = lgd.Units;
lgd.Units = 'normalized';
pos = lgd.Position;
center = pos(1:2) + pos(3:4) ./ 2;
pos(3:4) = pos(3:4) .* scale;
pos(1:2) = center - pos(3:4) ./ 2;
pos(1) = min(max(pos(1), 0), max(0, 1 - pos(3)));
pos(2) = min(max(pos(2), 0), max(0, 1 - pos(4)));
lgd.Position = pos;
lgd.Units = old_units;
end
