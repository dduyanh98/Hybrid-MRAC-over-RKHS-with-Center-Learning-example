function Box = box_from_centers(Xi, root_id)
Box.theta_L = min(Xi(:,1));
Box.theta_R = max(Xi(:,1));
Box.d_D = min(Xi(:,2));
Box.d_U = max(Xi(:,2));
Box.center = [(Box.theta_L + Box.theta_R)/2, (Box.d_D + Box.d_U)/2];
Box.corner_idx = [];
Box.depth = 0;
Box.active = true;
Box.root_id = root_id;
end
