function [Boxes, Xi_lib, y_lib] = refine_box(Boxes, boxID, Xi_lib, y_lib, p)
% REFINE_BOX
%   Split the leaf box Boxes(boxID) into 4 subboxes (depth+1).
%   Add midpoints and center to Xi_lib / y_lib if needed.
%   Parent box becomes inactive, children become active leaves.

B = Boxes(boxID);
if isfield(B, 'root_id')
    root_id = B.root_id;
else
    root_id = boxID;
end

% Do nothing if already at max depth
if B.depth >= p.max_depth
    return;
end

% Extract existing corner vertices in a fixed order:
% 1: BL, 2: BR, 3: TL, 4: TR
corner_idx = B.corner_idx;
C = Xi_lib(corner_idx, :);   % 4 x 2

BL = C(1,:);  BR = C(2,:);
TL = C(3,:);  TR = C(4,:);

% Compute midpoints in state-space
M_bottom = (BL + BR)/2;
M_left   = (BL + TL)/2;
M_right  = (BR + TR)/2;
M_top    = (TL + TR)/2;
M_center = (BL + TR)/2;

NewPoints = [M_bottom; M_left; M_right; M_top; M_center];
NewIdx    = zeros(5,1);

for k = 1:5
    pos = NewPoints(k,:);
    % Check if already exists (numerically)
    idx_exist = find( abs(Xi_lib(:,1)-pos(1)) < 1e-12 & ...
                      abs(Xi_lib(:,2)-pos(2)) < 1e-12, 1 );
    if ~isempty(idx_exist)
        NewIdx(k) = idx_exist;
    else
        Xi_lib = [Xi_lib; pos];
        try
            val = regressor_vector(pos');
        catch
            val = sin(pos(1));
        end
        y_lib = [y_lib; val];
        NewIdx(k) = size(Xi_lib,1);
    end
end

idx_MB = NewIdx(1);
idx_ML = NewIdx(2);
idx_MR = NewIdx(3);
idx_MT = NewIdx(4);
idx_MC = NewIdx(5);

% Easier access to scalar bounds
theta_L = B.theta_L;  theta_R = B.theta_R;
d_D     = B.d_D;      d_U     = B.d_U;

theta_mid = (theta_L + theta_R)/2;
d_mid     = (d_D     + d_U    )/2;

% Build 4 new child boxes (depth+1)
depth_child = B.depth + 1;
children = [];

% Child 1: bottom-left
b1.theta_L    = theta_L;
b1.theta_R    = theta_mid;
b1.d_D        = d_D;
b1.d_U        = d_mid;
b1.center     = [(theta_L+theta_mid)/2, (d_D+d_mid)/2];
b1.corner_idx = [corner_idx(1); idx_MB; idx_ML; idx_MC]; % BL, MB, ML, MC
b1.depth      = depth_child;
b1.active     = true;
b1.root_id    = root_id;
children = [children; b1];

% Child 2: bottom-right
b2.theta_L    = theta_mid;
b2.theta_R    = theta_R;
b2.d_D        = d_D;
b2.d_U        = d_mid;
b2.center     = [(theta_mid+theta_R)/2, (d_D+d_mid)/2];
b2.corner_idx = [idx_MB; corner_idx(2); idx_MC; idx_MR]; % MB, BR, MC, MR
b2.depth      = depth_child;
b2.active     = true;
b2.root_id    = root_id;
children = [children; b2];

% Child 3: top-left
b3.theta_L    = theta_L;
b3.theta_R    = theta_mid;
b3.d_D        = d_mid;
b3.d_U        = d_U;
b3.center     = [(theta_L+theta_mid)/2, (d_mid+d_U)/2];
b3.corner_idx = [idx_ML; idx_MC; corner_idx(3); idx_MT]; % ML, MC, TL, MT
b3.depth      = depth_child;
b3.active     = true;
b3.root_id    = root_id;
children = [children; b3];

% Child 4: top-right
b4.theta_L    = theta_mid;
b4.theta_R    = theta_R;
b4.d_D        = d_mid;
b4.d_U        = d_U;
b4.center     = [(theta_mid+theta_R)/2, (d_mid+d_U)/2];
b4.corner_idx = [idx_MC; idx_MR; idx_MT; corner_idx(4)]; % MC, MR, MT, TR
b4.depth      = depth_child;
b4.active     = true;
b4.root_id    = root_id;
children = [children; b4];

% Mark parent inactive, append children
Boxes(boxID).active = false;
Boxes = [Boxes; children];

end
