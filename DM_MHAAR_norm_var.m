function [outputs] = DM_MHAAR_norm_var(y, p, A, n, M, N, DP_eps, DP_delta, theta0, sigma_q, stats_name)


% This function infers the mean parameter of the normal distribution given
% the mean of x^p's for p in (0, 1)

beta = DP_eps/(2*log(2/DP_delta));
alfa = DP_eps/2 ;

theta = theta0;

Thetas = zeros(1, M);

% pseudorandom variables
z = randn(1, n);

A_sens = A^p;

dec_vec = zeros(1, M);

for u = 1:M
    %if mod(m, 1000) == 0
        %disp(m)
    %end
    theta_prop = theta + sigma_q*randn;
    if theta_prop > 0
        Z = [z; randn(N-1, n)];

        %transformation
        X = Z.*sqrt(theta);
        X_prop = Z.*sqrt(theta_prop);        
        
        X = min(A,max(-A,X));
        S = abs(X).^p;

        X_prop = min(A,max(-A, X_prop));
        S_prop = abs(X_prop).^p;
        
        if strcmp(stats_name, 'median') == 1
            var_y = smooth_sensitivity_median(S, beta, 0, A_sens);
            U = median(S , 2);
            var_y_prop = smooth_sensitivity_median(S_prop, beta, 0, A_sens);
            U_prop = median(S_prop, 2);
        
        elseif strcmp(stats_name, 'max') == 1
            var_y = smooth_sensitivity_max(S, beta, 0, A_sens);
            U = max(S, [], 2);
            var_y_prop = smooth_sensitivity_max(S_prop, beta, 0, A_sens);
            U_prop = max(S_prop, [], 2);
        end
        
        log_w = lap_log_pdf(y,U, var_y./alfa);
        log_sum = log_sum_exp(log_w);

        log_w_prop = lap_log_pdf(y, U_prop, var_y_prop./alfa);
        log_sum_prop = log_sum_exp(log_w_prop);

        log_r = log_sum_prop - log_sum;
        decision = rand < exp(log_r);
        
        if decision == 1
            theta = theta_prop;
            w = exp(log_w_prop - log_sum_prop);
            k = randsample(1:N, 1, 'true', w);
            z = Z(k, :);
            dec_vec(u) = decision;
        else
            w = exp(log_w - log_sum);
            k = randsample(1:N, 1, 'true', w);
            z = Z(k, :);
        end
    end
    Thetas(u) = theta;
end

outputs.Thetas = Thetas;
outputs.dec_vec = dec_vec;
