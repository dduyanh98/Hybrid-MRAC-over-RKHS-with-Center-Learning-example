function V = lyapunov_eps_value(y, p)
x = y(1:2);
x_ref = y(3:4);
e_tran = y(9:10);
eps = x - x_ref - e_tran;
V = eps' * p.P * eps;
end
