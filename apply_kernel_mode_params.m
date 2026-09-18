function p = apply_kernel_mode_params(p)
field_name = sprintf('kernel_mode%d_params', p.kernel_mode);

if ~isfield(p, field_name) || isempty(p.(field_name))
    return;
end

mode_params = p.(field_name);
param_fields = fieldnames(mode_params);
for i = 1:numel(param_fields)
    field = param_fields{i};
    p.(field) = mode_params.(field);
end
end
