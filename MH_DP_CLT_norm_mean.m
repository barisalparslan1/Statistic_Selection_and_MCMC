function [outputs] = MH_DP_CLT_norm_mean(y, theta0, a, n, A, DP_eps, M, sigma_q)

theta = theta0;
Delta = A.^a/n;

mu = theta*(a == 1) + (theta^3+ 3*theta)*(a == 3);
Sigma = (a == 1) + (9*theta^4 + 36*theta^2 + 15)*(a == 3);
U = Sigma/n + Delta^2/DP_eps;

theta = theta0;

Thetas = zeros(1, M);

for m = 1:M
    theta_prop = theta + sigma_q*randn;

    % Generate the sufficient statistics
    mu_prop = theta_prop*(a == 1) + (theta_prop^3 + 3*theta_prop)*(a == 3);
    Sigma_prop = (a == 1) + (9*theta_prop^4 + 36*theta_prop^2 + 15)*(a == 3);
    U_prop = Sigma_prop/n + Delta^2/DP_eps^2;

    % acceptance ratio
    log_p_prop = -0.5*log(2*pi*U_prop)-0.5*(y - mu_prop).^2/U_prop;
    log_p      = -0.5*log(2*pi*U)-0.5*(y - mu).^2/U;

    log_r = log_p_prop - log_p;
    
    decision = rand < exp(log_r);

    if decision == 1
        theta = theta_prop;

        mu = mu_prop;
        U = U_prop;
    end
    Thetas(m) = theta;
end

outputs.Thetas = Thetas;