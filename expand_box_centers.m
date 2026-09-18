function [Xi_sel, y_sel] = expand_box_centers(x_rep, Xi_current, Xi_lib, y_lib)
% EXPAND_BOX_CENTERS
%   Select four centers from a box twice as large as the current 4-corner box.
%   The current corner closest to x_rep is kept as the anchor corner, and the
%   expanded box extends away from that corner.

theta_vals = unique(Xi_current(:,1));
d_vals     = unique(Xi_current(:,2));

if numel(theta_vals) < 2 || numel(d_vals) < 2
    Xi_sel = Xi_current;
    y_sel = zeros(size(Xi_current,1),1);
    for i = 1:size(Xi_current,1)
        y_sel(i) = regressor_vector(Xi_current(i,:)');
    end
    return;
end

theta_L = min(theta_vals);
theta_R = max(theta_vals);
d_D     = min(d_vals);
d_U     = max(d_vals);

corners = [theta_L, d_D;
           theta_R, d_D;
           theta_L, d_U;
           theta_R, d_U];
[~, idx_near] = min(sum((corners - x_rep(:)').^2, 2));

w = theta_R - theta_L;
h = d_U - d_D;

switch idx_near
    case 1
        theta_new = [theta_L, theta_L + 2*w];
        d_new     = [d_D,     d_D     + 2*h];
    case 2
        theta_new = [theta_R - 2*w, theta_R];
        d_new     = [d_D,         d_D + 2*h];
    case 3
        theta_new = [theta_L,     theta_L + 2*w];
        d_new     = [d_U - 2*h, d_U];
    otherwise
        theta_new = [theta_R - 2*w, theta_R];
        d_new     = [d_U - 2*h, d_U];
end

Xi_target = [theta_new(1), d_new(1);
             theta_new(2), d_new(1);
             theta_new(1), d_new(2);
             theta_new(2), d_new(2)];

Xi_sel = zeros(4,2);
y_sel = zeros(4,1);
for i = 1:4
    pos = Xi_target(i,:);
    idx = find(abs(Xi_lib(:,1)-pos(1)) < 1e-12 & ...
               abs(Xi_lib(:,2)-pos(2)) < 1e-12, 1);
    if isempty(idx)
        [~, idx] = min(sum((Xi_lib - pos).^2, 2));
    end
    Xi_sel(i,:) = Xi_lib(idx,:);
    y_sel(i) = y_lib(idx);
end
end
