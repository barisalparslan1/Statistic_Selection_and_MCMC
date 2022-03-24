function [FIM] = FIM_DP_norm_var_CLT(theta, n, eps_DP, p, A)

% [FIM] = FIM_DP_norm_var_CLT(theta, n, eps_DP, p, A)
% 
% This function calculates the FIM of the noisy statistic with respect to the
% variance parameter of the normal distribution with mean 0.
% Inputs are 
% DP level: eps_DP, 
% data size: n, 
% moment of |x|: p, 
% boundary of |x|: A
% the variance parameter: theta.
%
% Last update: 7 March 2022

mu_deriv = theta^(p/2-1)*(p/2)*2^(p/2)*gamma((p+1)/2)/sqrt(pi);
S = (2*theta)^p*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2)^2/pi);
U = S/n + A^(2*p)/(n^2*eps_DP^2);
S_deriv = 2*p*(2*theta)^(p-1)*(gamma((2*p+1)/2)/sqrt(pi) - gamma((p+1)/2)^2/pi);
U_deriv = S_deriv/n;

FIM = mu_deriv^2/U + 0.5*(U_deriv^2)/(U^2);