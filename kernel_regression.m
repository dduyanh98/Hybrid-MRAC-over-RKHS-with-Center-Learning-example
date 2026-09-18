function [alpha, K_block, y_bar] = kernel_regression(f_Xi, Xi, y_data, sigma, lambda)
% KERNEL_REGRESSION  Build regression kernel and solve for alpha.
% Inputs:
%   f_Xi   : [N x d] feature points associated with centers (Soonyong style)
%   Xi     : [N x d] kernel centers
%   y_data : [N x m] outputs at centers (m=1 here)
%   sigma  : kernel bandwidth
%   lambda : nonnegative ridge regularization (scalar)
%
% Outputs:
%   alpha  : [(m*N) x 1]
%   K_block: [(m*N) x (m*N)]
%   y_bar  : [(m*N) x 1]

[N, d] = size(Xi); %#ok<NASGU>
m = size(y_data, 2);

K_block = zeros(m*N, m*N);
y_bar   = zeros(m*N, 1);

for i = 1:N
    for j = 1:N
        f_Xi_i = f_Xi(i,:);
        Xi_j   = Xi(j,:);
        k_val  = exp(-norm(f_Xi_i - Xi_j, 2)^2 / (2 * sigma^2));   % Gaussian

        for row_idx = 1:m
            for col_idx = 1:m
                row = m * (i-1) + row_idx;
                col = m * (j-1) + col_idx;
                if row_idx == col_idx
                    K_block(row, col) = k_val;
                else
                    K_block(row, col) = 0;
                end
            end
        end
    end

    % stack outputs
    for row_idx = 1:m
        row = m*(i-1) + row_idx;
        y_bar(row,1) = y_data(i, row_idx);
    end
end

% ridge (lambda * I)
if lambda > 0
    K_block = K_block + lambda * eye(size(K_block));
end

alpha = K_block \ y_bar;
end
