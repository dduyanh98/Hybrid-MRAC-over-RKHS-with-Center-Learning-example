function [UfgN, K_aprx, alpha] = kernel_prediction(x_state, x_centers, alpha, sigma, m)
% Predict function value at x_query using learned kernel model
% Inputs:
%   x_state     : [1 x d] input state to evaluate
%   x_centers   : [N x d] kernel centers
%   alpha       : [(m*N) x 1] learned coefficient vector
%   sigma       : kernel bandwidth
%   m           : output dimension (use m=1 here)
%
% Output:
%   UfgN        : [1 x m] predicted output vector

    [M, d] = size(x_state); %#ok<ASGLU>
    [N, ~] = size(x_centers);

    K_aprx = zeros(m, N*m);

    for i = 1:M         % M = 1
        for j = 1:N
            x_i  = x_state(i,:);
            Xi_i = x_centers(j,:);
            k_val = exp(-norm(x_i - Xi_i, 2)^2 / (2 * sigma^2));    % Gaussian Kernel

            for row_idx = 1:m
                for col_idx = 1:m
                    row = m*(i-1) + row_idx;
                    col = m*(j-1) + col_idx;
                    if row_idx == col_idx
                        K_aprx(row, col) = k_val;
                    end
                end
            end
        end
    end

    UfgN = K_aprx * alpha;
end
