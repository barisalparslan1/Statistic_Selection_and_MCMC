function [Thetas] = MHAAR_RB_DP_norm_var(y, theta0, DP_eps, A, p, sigma_q, M, N)

% [Thetas] = MHAAR_RB_DP_norm_mean(y, theta0, DP_eps, A, p, sigma_q, M, N)
% 
% This function is the MHAAR-RB algorithm for estimating the unknown var
% of a normal population with zero mean when the observations are
% released in privacy preserving noise.
% 
% Last update: 7 March 2022

n = length(y);
Delta = A^p;

lap_par_DP = Delta/DP_eps;

theta = theta0;
z = randn(1, n);

Thetas = zeros(1, M);

for m = 1:M
    % propose theta and the auxiliary variables
    theta_prop = theta + sigma_q*randn;

    if theta_prop > 0
        Z = [z; randn(N-1, n)];

        % calculate the weights
        X_par_0 = Z*sqrt(theta);
        log_w_0 = lap_log_pdf(y, abs(X_par_0).^p, lap_par_DP);
        w_max_0 = max(log_w_0, [], 1);
        log_sum_0 = log(sum(exp(log_w_0 - w_max_0))) + w_max_0;

        X_par_1 = Z*sqrt(theta_prop);
        log_w_1 = lap_log_pdf(y, abs(X_par_1).^p, lap_par_DP);
        w_max_1 = max(log_w_1, [], 1);
        log_sum_1 = log(sum(exp(log_w_1 - w_max_1))) + w_max_1;

        % calculate the acceptance ratio
        log_r = sum(log_sum_1 - log_sum_0);

        % decision and sampling of z
        if rand < exp(log_r)
            theta = theta_prop;
            W_temp = exp(log_w_1 - log_sum_1);
            temp_ind_vec = sum(rand(1, n) > cumsum(W_temp)) + 1;
            temp_mtx_ind_vec = N*(0:n-1) + temp_ind_vec;
            z = Z(temp_mtx_ind_vec);
        else
            W_temp = exp(log_sum_0 - log_sum_0);
            temp_ind_vec = sum(rand(1, n) > cumsum(W_temp)) + 1;
            temp_mtx_ind_vec = N*(0:n-1) + temp_ind_vec;
            z = Z(temp_mtx_ind_vec);
        end
    end

    % store the sample
    Thetas(m) = theta;

end
