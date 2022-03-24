function [outputs] = MH_DP_CLT_norm_var(y, theta0, p, n, A, DP_eps, M, sigma_q)

acceptance_counter = 0;
theta = theta0;

mu = theta.^(p/2).*2.^(p/2).*gamma((p+1)/2)/sqrt(pi);
S = (2*theta).^p.*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2).^2/pi);
U = S/n + A.^(2*p)/(n^2*DP_eps^2);

theta = theta0;

Thetas = zeros(1, M);
dec_vec = zeros(1,M);

for m = 1:M
    theta_prop = exp(log(theta)+ sigma_q*randn);

    % Generate the sufficient statistics
    mu_prop = theta_prop.^(p/2).*2.^(p/2).*gamma((p+1)/2)/sqrt(pi);
    S_prop = (2*theta_prop).^p.*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2).^2/pi);
    U_prop = S_prop/n + A.^(2*p)/(n^2*DP_eps^2);

    % acceptance ratio
    log_p_prop = -0.5*log(2*pi*U_prop)-0.5*(y - mu_prop).^2/U_prop;
    log_p      = -0.5*log(2*pi*U)-0.5*(y - mu).^2/U;

    log_r = log_p_prop - log_p;

    decision = rand < exp(log_r);

    if decision == 1
        theta = theta_prop;

        mu = mu_prop;
        U = U_prop;
        dec_vec(m) = decision;
    end
    Thetas(m) = theta;
end

outputs.Thetas = Thetas;
outputs.dec_vec = dec_vec;
