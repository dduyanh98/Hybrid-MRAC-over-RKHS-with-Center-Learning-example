function [Xi_sel, y_sel, BoxID] = select_centers_from_library(x_rep, Xi_lib, y_lib, Boxes, p)
% SELECT_CENTERS_FROM_LIBRARY
%   Given current state x_rep = [theta, theta_dot] (rad),
%   find the active leaf box containing it (or nearest) and
%   return its 4 corner centers as Xi_sel, y_sel.

theta = x_rep(1);
thetadot = x_rep(2);

% 1) Find all active boxes that contain the point
candidates = [];
for k = 1:numel(Boxes)
    if ~Boxes(k).active
        continue;
    end
    if theta >= Boxes(k).theta_L && theta <= Boxes(k).theta_R && ...
       thetadot >= Boxes(k).d_D && thetadot <= Boxes(k).d_U
        candidates = [candidates; k];
    end
end

if ~isempty(candidates)
    % If multiple (on boundaries), pick the deepest one
    depths = arrayfun(@(idx) Boxes(idx).depth, candidates);
    [~, idx_max] = max(depths);
    BoxID = candidates(idx_max);
else
    % 2) If no box contains the point (outside grid), pick nearest box center
    best_id = -1;
    best_d2 = inf;
    for k = 1:numel(Boxes)
        if ~Boxes(k).active
            continue;
        end
        c = Boxes(k).center;
        d2 = (theta - c(1))^2 + p.rho*(thetadot - c(2))^2;
        if d2 < best_d2
            best_d2 = d2;
            best_id = k;
        end
    end
    BoxID = best_id;
end

if BoxID < 0
    error('No active box found in select_centers_from_library.');
end

% 3) Use the four corner vertices of that box as centers
idx_corners = Boxes(BoxID).corner_idx;   % 4x1
Xi_sel = Xi_lib(idx_corners, :);
y_sel  = y_lib(idx_corners);

end
