function [mu,U] = mean_var_calc(theta,p,n)
    mu = theta.^(p/2).*2.^(p/2).*gamma((p+1)/2)/sqrt(pi);
    S = (2*theta).^p.*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2).^2/pi);
    U = S/n;
end