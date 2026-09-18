function W = lyapunov_tran_value(y, p)
e_tran = y(9:10);
W = e_tran' * p.P_tran * e_tran;
end
