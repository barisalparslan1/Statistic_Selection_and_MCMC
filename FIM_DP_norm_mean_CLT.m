function [FIM] = FIM_DP_norm_mean_CLT(theta, n, eps_DP, p, A)

% [FIM] = FIM_DP_norm_var_CLT(theta, n, eps_DP, p, A)
%
% This function calculates the FIM of the noisy statistic with respect to the
% mean parameter of the normal distribution with variance 1.
% Inputs are
% DP level: eps_DP,
% data size: n,
% moment of |x|: p,
% boundary of x: (0, A)
% the mean parameter: theta.
%
% Last update: 7 March 2022

if p == 1
    mu_deriv = 1;
    U = 1/n + A.^2/(n^2*eps_DP^2);
    U_deriv = 0;
    FIM = mu_deriv.^2./U + 0.5*(U_deriv.^2)./(U.^2);
elseif p == 3
    mu_deriv = 3*theta^2 + 3;
    S = 9*theta^4 + 36*theta^2 + 15;
    U = S/n + A.^(2*3)/(n^2*eps_DP^2);
    U_deriv = (36*theta^3 + 72*theta)/n;
    FIM =  mu_deriv.^2./U + 0.5*(U_deriv.^2)./(U.^2);
end