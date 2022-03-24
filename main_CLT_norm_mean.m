% FIM calculation for Example 1 in the manuscript
%
% These are the experiments regarding additive statistics. Gaussian
% mechanism is used to provide eps, delta -DP.
%
% The population is normal with zero mean and unknonw variance theta.
% We compare x^p for different p as statistics in terms of FIM and MSE.
%
% The MCMC algorithm is a MH algorithm.
%
% Use the first code block for plots in Example 1 in the manuscript
%
% Last update: 22 March 2022

close all; clc; clear;
%% FIM values in Example 1

A = 10;
n = 100; p_vec = [1 3]; L_p = length(p_vec);
eps_DP_vec = [1 inf]; L_eps = length(eps_DP_vec);

Theta = 0.1:0.1:10; L_theta = length(Theta);

FIM = zeros(L_eps, L_p, L_theta);

for i1 = 1:L_eps
    eps_DP = eps_DP_vec(i1);
    for i2 = 1:L_p
        p = p_vec(i2);
        for i3 = 1:L_theta
            theta = Theta(i3);
            FIM(i1, i2, i3) = FIM_DP_norm_mean_CLT(theta, n, eps_DP, p, A);
        end
    end
end

for i = 1:L_eps
    subplot(1, L_eps, i);
    plot(Theta, log(squeeze(FIM(i, 1, :))), 'b', Theta, log(squeeze(FIM(i, 2, :))), 'r');
    set(gcf, 'Position',  [100, 100, 200, 150]);
    xlabel('$\theta$', 'Interpreter', 'latex');
    ylabel('$\log F(\theta)$', 'Interpreter', 'latex');
    title(['$\epsilon$' sprintf( '= %.2f', eps_DP_vec(i))], 'Interpreter', 'Latex');
    filenametoprint = sprintf('Normal_mean_A_%d_n_%d_eps_%02d', A, n, 100*eps_DP_vec(i));
    legend('$a = 1$', '$a = 3$', 'Interpreter', 'latex', 'Location','southeast');
    print(gcf,'-depsc2', filenametoprint);   
end


%% This is the main code for trying the ABC_DP algorithms

clear; clc; close all; fc = 0;

theta = 9.5;

A = 10;
n = 100;

DP_eps = 0.1;
X = randn(1,n) + theta;


M = 100000;
theta0 = 1;

sigma_q = sqrt(1/n);

a_vec = [1 3]; L_a = length(a_vec);
outputs = cell(1,L_a);

tic;
for i = 1:L_a
    a = a_vec(i);

    Delta = A^a/n;
    v = randn*sqrt(Delta^2/DP_eps^2);
    y = mean(X.^a) + v;

    [outputs{i}] = MH_DP_CLT_norm_mean(y, theta0, a, n, A, DP_eps, M, sigma_q);
end
toc;

fc = fc + 1; figure(fc);
subplot(2, 1, 1);
plot(outputs{1}.Thetas); hold on; plot(outputs{2}.Thetas); hold off; legend('a = 1','a = 3');

subplot(2, 2, 3);
histogram(outputs{1}.Thetas(end/4:end), 50);
subplot(2, 2, 4);
histogram(outputs{2}.Thetas(end/4:end), 50);

mean(outputs{1}.Thetas(end/4:end).^2);
mean(outputs{2}.Thetas(end/4:end).^2);

