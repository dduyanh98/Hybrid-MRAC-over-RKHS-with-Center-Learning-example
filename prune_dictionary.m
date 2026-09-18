function [Xi_lib, y_lib, BoxCenters, BoxCornerIdx] = prune_dictionary( ...
    x_rep, Xi_lib, y_lib, BoxCenters, BoxCornerIdx, t_event, p)
% Removes boxes whose centers are TOO FAR from current state.
%
% Inputs:
%   x_rep         : [1x2] current state (rad)
%   Xi_lib        : Nx2 remaining dictionary points
%   y_lib         : Nx1 regression outputs
%   BoxCenters    : [B x 2] box centers (rad)
%   BoxCornerIdx  : cell array {B}, each cell = [4x1] indices of Xi_lib
%   t_event       : event time (sec)
%   p             : parameter struct with fields:
%                       p.r0    -- initial prune radius
%                       p.gamma -- growth rate (per sec)
%                       p.rho   -- velocity weight (recommended 0.25)

theta = x_rep(1);
thetadot = x_rep(2);

% dynamic threshold
rmax = p.r0 + p.gamma * t_event;

% mask of boxes to KEEP
keepBox = true(length(BoxCenters),1);

for b = 1:length(BoxCenters)
    theta_c    = BoxCenters(b,1);
    thetadot_c = BoxCenters(b,2);

    % weighted distance
    d2 = (theta - theta_c)^2 + p.rho * (thetadot - thetadot_c)^2;
    d  = sqrt(d2);

    if d > rmax
        keepBox(b) = false;  % mark box for removal
    end
end

% Indices of boxes to remove
removeBoxes = find(~keepBox);

% Gather all dictionary indices to remove
rm_idx = [];
for rb = removeBoxes'
    rm_idx = [rm_idx; BoxCornerIdx{rb}(:)];
end
rm_idx = unique(rm_idx);

% Now remove them from Xi_lib and y_lib
Xi_lib(rm_idx,:) = [];
y_lib(rm_idx)    = [];

% Remove those boxes
BoxCenters(removeBoxes,:) = [];
BoxCornerIdx(removeBoxes) = [];

end
