function [outputs] = MHAAR_DP_CLT_Lap_norm_var(y, theta0,mu0, p, n, A, DP_eps, M, sigma_q,N)

theta = theta0;
u = mu0;
sigma = A.^(p)/(n*DP_eps);
Thetas = zeros(1, M);

for m = 1:M
    
    theta_prop = exp(log(theta)+ sigma_q*randn);
    theta_for_q = max(theta_prop,theta);
    
    [mu,U] = norm_variance_mean_var(theta,p,n);
    [mu_for_q,U_for_q] = norm_variance_mean_var(theta_for_q,p,n);
    [mu_prop,U_prop] = norm_variance_mean_var(theta_prop,p,n);
    
    u_js = mu_for_q+sqrt(U_for_q)*randn(1, N);
    u_js(1)=u;
       
    % Generate the sufficient statistics --- for theta
    log_w = lap_log_pdf(y,u_js,sigma)-0.5.*(u_js - mu).^2./U ...
        - 0.5*log(2*pi*U) ...
        + 0.5.*(u_js - mu_for_q).^2./U_for_q...
        + 0.5*log(2*pi*U_for_q);
    log_w_sum = log_sum_exp(log_w);
    log_Z_T = log_w_sum - log(N);

    % Generate the sufficient statistics --- for theta prime
    
    log_w_prop = lap_log_pdf(y,u_js,sigma)-0.5.*(u_js - mu_prop).^2./U_prop...
        -0.5*log(2*pi*U_prop) + 0.5*(u_js - mu_for_q).^2./U_for_q + 0.5*log(2*pi*U_for_q);
    log_w_sum_prop = log_sum_exp(log_w_prop);
    log_Z_T_prime = log_w_sum_prop - log(N);
    
    log_accept_rate = log_Z_T_prime - log_Z_T;
    decision = rand < exp(log_accept_rate);

    if decision == 1
        w_prop = exp(log_w_prop - log_w_sum_prop);
        k = randsample(1:N, 1, true, w_prop);
        u = u_js(k);
        theta = theta_prop;
    else
        w = exp(log_w - log_w_sum);
        k = randsample(1:N, 1, true, w);
        u = u_js(k);
    end
    Thetas(m) = theta;
end
outputs.Thetas = Thetas;
end