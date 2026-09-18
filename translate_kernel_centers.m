function [Xi_new, y_new] = translate_kernel_centers(Xi_current, old_anchor, new_anchor)
% Shift the active kernel stencil without changing its size or center count.
if isempty(Xi_current) || isempty(old_anchor) || isempty(new_anchor)
    Xi_new = Xi_current;
else
    delta = new_anchor(:)' - old_anchor(:)';
    Xi_new = Xi_current + delta;
end

y_new = zeros(size(Xi_new,1), 1);
for i = 1:size(Xi_new,1)
    y_new(i) = regressor_vector(Xi_new(i,:)');
end
end
