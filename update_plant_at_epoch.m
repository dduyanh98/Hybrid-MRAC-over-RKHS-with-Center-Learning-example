function p = update_plant_at_epoch(p)
% UPDATE_PLANT_AT_EPOCH  -- change physical pendulum parameters (m,l,k)
% at each epoch jump and recompute reference model & Lyapunov matrix.
%
% This version:
%   - keeps kp, kd as SCALARS
%   - uses a simple deterministic schedule for (m,l,k)
%   - recomputes A_ref and P cleanly every epoch

    % ------------------------------------------------------------------
    % Epoch index bookkeeping
    % ------------------------------------------------------------------
    if ~isfield(p,'epoch_index') || isempty(p.epoch_index)
        p.epoch_index = 0;
    end
    p.epoch_index = p.epoch_index + 1;

    if ~isfield(p,'plant_schedule') || isempty(p.plant_schedule)
        plant_schedule = 'current';
    else
        plant_schedule = p.plant_schedule;
    end

    switch plant_schedule
        case 'current'
            idx = mod(p.epoch_index-1, 5) + 1;

            switch idx
                case 1
                    p.m = 1.40;  p.l = 0.75;   p.k = 0.55;
                case 2
                    p.m = 1.15;  p.l = 0.95;   p.k = 0.42;
                case 3
                    p.m = 1.55;  p.l = 0.70;   p.k = 0.65;
                case 4
                    p.m = 1.25;  p.l = 0.88;   p.k = 0.48;
                case 5
                    p.m = 1.35;  p.l = 0.82;   p.k = 0.58;
            end

        case 'aggressive'
            idx = mod(p.epoch_index-1, 5) + 1;

            switch idx
                case 1
                    p.m = 7.3;   p.l = 0.08;   p.k = 0.05;
                case 2
                    p.m = 0.06;  p.l = 11.2;   p.k = 0.1;
                case 3
                    p.m = 9.8;   p.l = 0.06;   p.k = 8.0;
                case 4
                    p.m = 2.6;   p.l = 5.3;    p.k = 7.2;
                case 5
                    p.m = 0.15;  p.l = 40;     p.k = 13.9;
            end

        otherwise
            error('Unknown p.plant_schedule = %s. Use ''current'' or ''aggressive''.', plant_schedule);
    end

    % ------------------------------------------------------------------
    % Recompute reference-model gains based on new length l
    % ------------------------------------------------------------------
    if ~isfield(p,'g') || isempty(p.g)
        p.g = 9.81;
    end

    wn   = sqrt(p.g / p.l_est);
zeta = 0.75;

    p.kp = wn^2;              % SCALAR
    p.kd = 2*zeta*wn;         % SCALAR

    % ------------------------------------------------------------------
    % Recompute transient dynamics and Lyapunov matrices
    % ------------------------------------------------------------------
    if ~isfield(p,'Q') || isempty(p.Q)
        p.Q = eye(2);
    end

    A_ref = [0 1;
            -p.kp -p.kd];     % 2x2, no dimension mismatch

    if ~isfield(p,'Q_tran') || isempty(p.Q_tran)
        p.Q_tran = eye(2);
    end

    slowest_ref = max(real(eig(A_ref)));
    p.A_tran = diag([2.0, 2.5] * slowest_ref);
    p.P = lyap(p.A_tran', p.Q);
    p.P_tran = lyap(p.A_tran', p.Q_tran);

end
