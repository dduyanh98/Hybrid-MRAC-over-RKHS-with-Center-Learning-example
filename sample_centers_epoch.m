function [Xi_new, f_Xi_new, y_data_new] = sample_centers_epoch(x_traj, p, Xi_prev)
% SAMPLE_CENTERS_EPOCH - pick one new center from recent trajectory samples.
%   x_traj  : [K x d] trajectory samples, here [theta, theta_dot]
%   p       : parameter struct, uses p.r_min as the separation radius
%   Xi_prev : [N x d] existing centers, may be empty
%
% The candidate is the trajectory sample whose nearest existing center is
% farthest away. It is accepted only if that distance is larger than p.r_min.

assert(size(x_traj,2)==2,'Expecting 2D state [theta, theta_dot].');

if isempty(x_traj)
    Xi_new = [];
    f_Xi_new = [];
    y_data_new = [];
    return;
end

if isempty(Xi_prev)
    Xi_new = x_traj(end,:);
else
    dists_to_centers = zeros(size(x_traj,1), size(Xi_prev,1));
    for j = 1:size(Xi_prev,1)
        dists_to_centers(:,j) = sqrt(sum((x_traj - Xi_prev(j,:)).^2, 2));
    end

    nearest_center_dist = min(dists_to_centers, [], 2);
    [d_max, idx_max] = max(nearest_center_dist);

    if d_max > p.r_min
        Xi_new = x_traj(idx_max,:);
    else
        Xi_new = [];
    end
end

if isempty(Xi_new)
    f_Xi_new = [];
    y_data_new = [];
    return;
end

f_Xi_new = Xi_new;
try
    y_data_new = regressor_vector(Xi_new');
catch
    y_data_new = sin(Xi_new(1));
end
end
