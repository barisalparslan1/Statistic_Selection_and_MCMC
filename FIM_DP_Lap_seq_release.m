function FIM = FIM_DP_Lap_seq_release(theta, p,n, A, DP_eps, N, M)

log_deriv_p = zeros(1,N);

sigma = A.^(p)/(DP_eps);
for i = 1:N
    x_for_y = sqrt(theta)*randn;
    v = laprnd(1,1,0,sigma);
    y = abs(x_for_y^p) + v;
    
    X = sqrt(theta)*randn(1,M);
    log_w = lap_log_pdf(y, abs(X.^p), sigma); 
    
    G_sub = -0.5/(theta) + 0.5.*X.^2/(theta^2);
    w = exp(log_w - log_sum_exp(log_w));
    
    log_deriv_p(i) = sum(G_sub.*w);
end
FIM  = n*sum(log_deriv_p.^2)/N;
end