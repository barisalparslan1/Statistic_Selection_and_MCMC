function FIM = FIM_DP_CLT_Lap_norm_var(theta, p, n, A, DP_eps, N, M)

log_deriv_p = zeros(1,M);

[mu,U] = norm_variance_mean_var(theta,p,n);

mu_deriv = theta.^(p/2-1).*(p/2).*2.^(p/2).*gamma((p+1)/2)/sqrt(pi);
S_deriv = 2*p.*(2*theta).^(p-1).*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2).^2/pi);
U_deriv = S_deriv/n;


sigma = A.^(p)/(n*DP_eps);
for i = 1:M
    y = mu+sqrt(U)*randn + laprnd(1,1,0,sigma); 
    u_js = mu+sqrt(U)*randn(1, N);
    log_w = lap_log_pdf(y,u_js,sigma);
    G = -1/2.*(U.*2.*pi).^(-3/2)*U_deriv + 1/2.*(u_js-mu).^2/U.^2*U_deriv +(u_js-mu).*mu_deriv./U;
    w = exp(log_w - log_sum_exp(log_w));
    log_deriv_p(i) = sum(G.*w);
    
end
FIM  = sum(log_deriv_p.^2)/M;
end

%G = -1/2*((2.*(u_js-mu).*(-(p/2).*theta.^(p/2-1)*2.^(p/2).*gamma((p+1)/2)/sqrt(pi))).*U...
    %-(p.*theta.^(p-1).*2.^p.*(gamma((2*p+1)/2)/sqrt(pi)...
    %- gamma((p+1)/2).^2/pi)).*(u_js-mu).^2)/(U.^2);
