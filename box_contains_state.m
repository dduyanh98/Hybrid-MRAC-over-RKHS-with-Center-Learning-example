function tf = box_contains_state(Box, x_rep)
tol = 1e-12;
theta = x_rep(1);
thetadot = x_rep(2);

tf = theta >= Box.theta_L - tol && theta <= Box.theta_R + tol && ...
     thetadot >= Box.d_D - tol && thetadot <= Box.d_U + tol;
end
