function [outputs] = DM_MHAAR_norm_mean(y, p, A, n, M, N, DP_eps, theta0, sigma_q)

% This function infers the mean parameter of the normal distribution given
% the mean of x^p's for p in (0, 1)

theta = theta0;
acceptance_counter = 0;

Thetas = zeros(1, M);
var_y = (A^p/DP_eps/n)^2;

% pseudorandom variables
z = randn(1, n);
for m = 1:M
    %if mod(m, 1000) == 0
        %disp(m)
    %end
    theta_prop = theta + randn*sigma_q;
    
    Z = [z; randn(N-1, n)];
    
    %transformation
    X = theta+Z;
    X_prop = theta_prop+Z;
    
    U = mean(abs(X).^p.*sign(X), 2);
    U_prop = mean(abs(X_prop).^p.*sign(X_prop), 2);
    
    log_w = -(y - U).^2/(2*var_y);
    log_w_prop = -(y - U_prop).^2/(2*var_y);

    log_sum = log_sum_exp(log_w);
    log_sum_prop = log_sum_exp(log_w_prop);
    
    log_r = log_sum_prop - log_sum;
    decision = rand < exp(log_r);
    if decision == 1
        theta = theta_prop;
        w = exp(log_w_prop - log_sum_prop);
        k = randsample(1:N, 1, 'true', w);
        z = Z(k, :);
        if m>=M/4  % count acceptance ratio after burn-in
            acceptance_counter = acceptance_counter + 1;
        end
    else
        w = exp(log_w - log_sum);
        k = randsample(1:N, 1, 'true', w);
        z = Z(k, :);
    end
    Thetas(m) = theta;
end
acceptance_ratio = 4*acceptance_counter / (3*M);

outputs.Thetas = Thetas;
outputs.acceptance_ratio = acceptance_ratio;
end

