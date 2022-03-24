function [outputs] = PMMH_DP_CLT_Lap_norm_var(y, theta0, p, n, A, DP_eps, M, sigma_q, N)

theta = theta0;

[mu,U] = norm_variance_mean_var(theta,p,n);

sigma = A.^(p)/(n*DP_eps);

u_js = mu+sqrt(U)*randn(1, N);
log_w = lap_log_pdf(y,u_js,sigma);
log_w_sum = log_sum_exp(log_w);
log_Z_T = log_w_sum - log(N);

Thetas = zeros(1, M);
dec_vec = zeros(1, M);

for m = 1:M
    theta_prop = exp(log(theta)+ sigma_q*randn);

    % Generate the sufficient statistics
    [mu_prop,U_prop] = norm_variance_mean_var(theta_prop,p,n);

    u_js = mu_prop+sqrt(U_prop)*randn(1, N);
    log_w = lap_log_pdf(y,u_js,sigma);
    log_w_sum = log_sum_exp(log_w);
    log_Z_T_prime = log_w_sum - log(N);

    % acceptance probability
    log_accept_rate = log_Z_T_prime - log_Z_T;

    decision = rand < exp(log_accept_rate);

    if decision == 1
        theta = theta_prop;
        log_Z_T = log_Z_T_prime;
        dec_vec(m) = decision;
    end
    Thetas(m) = theta;
end

outputs.Thetas = Thetas;
outputs.dec_vec = dec_vec;


