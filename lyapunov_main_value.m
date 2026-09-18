function V = lyapunov_main_value(y, p)
x = y(1:2);
x_ref = y(3:4);
Kx = y(5:6);
Kr = y(7);
e_tran = y(9:10);
Kg = y(14:15);
Theta_rkhs = y(p.idx_theta_rkhs);

e = x - x_ref;
eps = e - e_tran;

V = eps' * p.P * eps ...
    + Kx' * (p.Gamma_x \ Kx) ...
    + Kr * (p.Gamma_r \ Kr);

if p.architecture_mode == 3
    V = V + Kg' * (p.Gamma_g \ Kg);
end

if p.architecture_mode >= 2
    V = V + (Theta_rkhs' * Theta_rkhs) / p.Gamma_f;
end
end
