function [Xi_sel, y_sel] = quadtree_grid_centers(Box, grid_n)
% Build a local 2x2, 4x4, or 6x6 grid by adding uniform outer layers
% around the active dictionary box.
grid_n = max(2, 2 * round(grid_n / 2));
theta_step = Box.theta_R - Box.theta_L;
d_step = Box.d_U - Box.d_D;
layers = (grid_n - 2) / 2;

theta_vals = Box.theta_L + (-layers:(layers+1)) * theta_step;
d_vals = Box.d_D + (-layers:(layers+1)) * d_step;

[TH, DD] = ndgrid(theta_vals, d_vals);
Xi_sel = [TH(:), DD(:)];

y_sel = zeros(size(Xi_sel,1), 1);
for i = 1:size(Xi_sel,1)
    y_sel(i) = regressor_vector(Xi_sel(i,:)');
end
end
