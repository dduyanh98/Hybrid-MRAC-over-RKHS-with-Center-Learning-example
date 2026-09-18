function [Xi_sel, y_sel] = full_box_grid_centers(Box, grid_n)
% Build all centers inside a box: 2x2, 3x3, or 5x5.
grid_n = max(2, round(grid_n));

theta_vals = linspace(Box.theta_L, Box.theta_R, grid_n);
d_vals = linspace(Box.d_D, Box.d_U, grid_n);

[TH, DD] = ndgrid(theta_vals, d_vals);
Xi_sel = [TH(:), DD(:)];

y_sel = zeros(size(Xi_sel,1), 1);
for i = 1:size(Xi_sel,1)
    y_sel(i) = regressor_vector(Xi_sel(i,:)');
end
end
