function [Xseg, p] = accumulate_positive_budget_segment(Xseg, p)
S_eps = Xseg(1,p.idx_S_eps);
S_tran = Xseg(1,p.idx_S_tran);

for k = 2:size(Xseg,1)
    X_prev = Xseg(k-1,:)';
    X_now = Xseg(k,:)';

    Veps_prev = lyapunov_eps_value(X_prev, p);
    Veps_now = lyapunov_eps_value(X_now, p);
    Vtran_prev = lyapunov_tran_value(X_prev, p);
    Vtran_now = lyapunov_tran_value(X_now, p);

    if Veps_now > Veps_prev
        S_eps = S_eps + (Veps_now - Veps_prev);
    end

    if Vtran_now > Vtran_prev
        S_tran = S_tran + (Vtran_now - Vtran_prev);
    end

    Xseg(k,p.idx_S_eps) = S_eps;
    Xseg(k,p.idx_S_tran) = S_tran;
end

p.ref_budget_cum = S_eps;
p.tran_budget_cum = S_tran;
end
