function FIM = FIM_DP_norm_var(theta, p, n, A, DP_eps, DP_delta, N, M, stats_name)

log_deriv_p = zeros(1,N);

beta = DP_eps/(2*log(2/DP_delta));
alfa = DP_eps/2 ;
A_sens = A^p;

for i = 1:N

    x_for_y = sqrt(theta)*randn(1,n);
    x_for_y = min(A,max(-A,x_for_y));
    S = abs(x_for_y).^p;
    
    if strcmp(stats_name, 'median') == 1
        sigma_y = smooth_sensitivity_median(S, beta, 0, A_sens);
        y = median(S) + laprnd(1,1,0,sigma_y/alfa);
    elseif strcmp(stats_name, 'max') == 1
        sigma_y = smooth_sensitivity_max(S, beta, 0, A_sens);
        y = max(S) + laprnd(1,1,0,sigma_y/alfa);
    end

    X = sqrt(theta)*randn(M,n);
    X = min(A,max(-A,X));
    S = abs(X).^p;

    if strcmp(stats_name, 'median') == 1
        sigma = smooth_sensitivity_median(S, beta, 0, A_sens);
        log_w = lap_log_pdf(y, median(S, 2), sigma/alfa);
    elseif strcmp(stats_name, 'max') == 1
        sigma = smooth_sensitivity_max(S, beta, 0, A_sens);
        log_w = lap_log_pdf(y, max(S, [], 2), sigma/alfa);
    end
    G_sub = sum(-0.5/(theta) + 0.5.*X.^2/(theta^2), 2);
    w = exp(log_w - log_sum_exp(log_w));

    log_deriv_p(i) = sum(G_sub.*w);    
end

FIM  = sum(log_deriv_p.^2)/N;