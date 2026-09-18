function s = find_s(s_previous, weighted_e_squared, alpha)

% This function allows to compute s such that \sum s is convergent; see the
% comments after (40) in the paper

% temporary_s = s+1;
% 
% while series_element(temporary_s) > weighted_e_squared
%     temporary_s = temporary_s + 1;
% end

if nargin < 3 || isempty(alpha)
    alpha = 1.001;
end
s = max(ceil(weighted_e_squared^(-1/alpha)),s_previous + 1);
