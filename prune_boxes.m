function Boxes = prune_boxes(x_rep, Boxes, t_event, p, ActiveBoxID)
% PRUNE_BOXES
%   Deactivate boxes (leaves) whose centers are too far from x_rep,
%   except the currently active box.

theta    = x_rep(1);
thetadot = x_rep(2);

rmax = p.r0 + p.gamma * t_event;
if rmax < 0
    rmax = 0;
end

for k = 1:numel(Boxes)
    if ~Boxes(k).active
        continue;
    end
    if k == ActiveBoxID
        continue;   % never prune the box we are using
    end

    c = Boxes(k).center;
    d2 = (theta - c(1))^2 + p.rho*(thetadot - c(2))^2;
    d  = sqrt(d2);

    if d > rmax
        Boxes(k).active = false;
    end
end
end
