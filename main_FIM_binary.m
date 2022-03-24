% FIM calculation for Example 4 in the manuscript
%
% The population distribution is Bernoulli with the parameter of theta. Gaussian mechanism
% and randomised binary responses are used to provide eps, delta -DP.
% 
% This is the experiment regarding binary responses. 
%
% Use the first code block for obtaining FIM values for Example 4.
%
% Last update: 22 March 2022

clc; clear; close all; fc = 0;
%% FIM values and plotting in the Example 4
n = 100;
DP_eps_vec = (0.01:0.01:50)';
theta_vec = 0.001:0.001:0.999;

alpha = (exp(DP_eps_vec)-1)./(exp(DP_eps_vec)+1);
tau = alpha*theta_vec + 1./(exp(DP_eps_vec)+1);
FIM1 = n*alpha.^2./(tau.*(1 - tau));

Var_Z = theta_vec.*(1 - theta_vec)/n + 1./(DP_eps_vec.^2)/n;
Var_Z_deriv = (1 - 2.*theta_vec)/n;
FIM2 = (Var_Z_deriv.^2)./(Var_Z.^2) + 1./Var_Z;

Var_Z = theta_vec.*(1 - theta_vec)/n + 1./(DP_eps_vec.^2)/n^2;
Var_Z_deriv = (1 - 2.*theta_vec)/n;
FIM3 = (Var_Z_deriv.^2)./(Var_Z.^2) + 1./Var_Z;

fc = fc + 1; figure(fc);
% subplot(2, 2, 1);
subplot('Position', [0.1, 0.6, 0.40, 0.35]);

imagesc(theta_vec, DP_eps_vec, log(1 + abs(FIM2 - FIM1)).*sign(FIM2 - FIM1)); colormap(gray); colorbar;
% imagesc(theta_vec, DP_eps_vec, FIM2 - FIM1); colormap(gray); colorbar;
title('$\log(1 + |F_{2}(\theta) - F_{1}(\theta)|)$ sign($F_{2}(\theta) - F_{1}(\theta)$)', 'Interpreter', 'Latex');
xlabel('$\theta$', 'Interpreter', 'Latex');
ylabel('$\epsilon$', 'Interpreter', 'Latex');

subplot('Position', [0.55, 0.6, 0.40, 0.35]);
imagesc(theta_vec, DP_eps_vec, FIM2 - FIM1 > 0); colormap(gray); % colorbar('XTickLabel',{'0','1'}, 'XTick', [0, 1]);
title('Black: $F_{2}(\theta) < F_{1}(\theta)$, White: $F_{2}(\theta) > F_{1}(\theta)$', 'Interpreter', 'Latex');
xlabel('$\theta$', 'Interpreter', 'Latex');
ylabel('$\epsilon$', 'Interpreter', 'Latex');

subplot('Position', [0.1, 0.1, 0.40, 0.35]);
imagesc(theta_vec, DP_eps_vec, log(1 + abs(FIM3 - FIM1)).*sign(FIM3 - FIM1)); colormap(gray); colorbar;
% imagesc(theta_vec, DP_eps_vec, FIM3 - FIM1); colormap(gray); colorbar;
title('$\log(1 + |F_{3}(\theta) - F_{1}(\theta)|)$ sign($F_{3}(\theta) - F_{1}(\theta)$)', 'Interpreter', 'Latex');
xlabel('$\theta$', 'Interpreter', 'Latex');
ylabel('$\epsilon$', 'Interpreter', 'Latex');

subplot('Position', [0.55, 0.1, 0.40, 0.35]);
imagesc(theta_vec, DP_eps_vec, FIM3 - FIM1 > 0); colormap(gray); % colorbar('XTickLabel',{'0','1'}, 'XTick', [0, 1]);
title('Black: $F_{3}(\theta) < F_{1}(\theta)$, White: $F_{3}(\theta) > F_{1}(\theta)$', 'Interpreter', 'Latex');
xlabel('$\theta$', 'Interpreter', 'Latex');
ylabel('$\epsilon$', 'Interpreter', 'Latex');


%%
n_vec = 10:100;
theta = 0.25;
DP_eps = 1;

alpha = (exp(DP_eps)-1)./(exp(DP_eps)+1);
tau = alpha*theta+ 1./(exp(DP_eps)+1);
FIM1 = n_vec.*alpha.^2./(tau.*(1 - tau));

Var_Z = theta.*(1 - theta)./n_vec + 1./(DP_eps.^2)./n_vec;
Var_Z_deriv = (1 - 2.*theta)./n_vec;
FIM2 = (Var_Z_deriv.^2)./(Var_Z.^2) + 1./Var_Z;

Var_Z = theta.*(1 - theta)./n_vec + 1./(DP_eps.^2)./(n_vec.^2);
Var_Z_deriv = (1 - 2.*theta)./n_vec;
FIM3 = (Var_Z_deriv.^2)./(Var_Z.^2) + 1./Var_Z;

fc = fc + 1; figure(fc);
plot(n_vec, log(FIM1), 'b', n_vec, log(FIM2), 'r', n_vec, log(FIM3), 'k');
legend({'$F_{1}(\theta)$', '$F_{2}(\theta)$', '$F_{3}(\theta)$'}, 'Interpreter','Latex');
xlabel('n');
ylabel('(log-) FIM');
title('Fisher information vs $n$ - $\epsilon = 1$, $\theta = 0.25$', 'Interpreter','Latex');
% set(gca, 'xticklabel', DP_eps_vec);
% set(gca, 'yticklabel', DP_eps_vec);