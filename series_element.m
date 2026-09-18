function element = series_element(s, alpha)

% element = 1/s^2;

% k = 3;
% element = 1/k^s;

if nargin < 2 || isempty(alpha)
    alpha = 1.001;
end
element = 1./s.^alpha;
