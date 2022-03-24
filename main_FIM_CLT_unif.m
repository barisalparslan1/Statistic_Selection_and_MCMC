% FIM calculation for Example 3 in the manuscript
%
% The population is uniform with width of 2*theta. Gaussian mechanism is
% used to provide eps, delta -DP.
%
% We compare |x|^p for different p as statistics in terms of FIM.
%
% Use the first code block for obtaining FIM values for Example 3
%
% Last update: 22 March 2022

close all; clc; clear;
%% FIM Values in Example 3

Theta = 0.01:0.01:2;
L_theta = length(Theta);
A = 10;
eps_DP_vec = [1 inf]; L_eps = length(eps_DP_vec);
a_vec = [0.01 0.1 0.5 1 2]; L_a = length(a_vec);
n = 100;
FIM = cell(1, L_eps);

for i = 1:L_eps
    eps_DP = eps_DP_vec(i);
    for j = 1:L_theta
        theta = Theta(j);
        mu_deriv = a_vec.*(theta).^(a_vec-1)./(a_vec+1);
        U = theta.^(2*a_vec).*(a_vec.^2./((a_vec+1).^2.*(2*a_vec + 1)))/n + A.^(2*a_vec)/(n^2*eps_DP^2);
        U_deriv = theta.^(2*a_vec-1).*(2*a_vec.^3./((a_vec+1).^2.*(2*a_vec + 1)))/n;
        temp = mu_deriv.^2./U + 0.5*(U_deriv.^2)./(U.^2);
        for k = 1:L_a
            FIM{k}(i, j) = temp(k);
        end
    end
end

legend_labels = cell(1, L_a);
for i = 1:L_eps
    subplot(1, L_eps, i)
    for k = 1:L_a
        a_temp = a_vec(k);
        legend_labels{k} = sprintf('$a$ = %.2f', a_temp);
        plot(Theta, log(FIM{k}(i, :)));
        hold on;
    end
    hold off;

    set(gcf, 'Position',  [100, 100, 320, 240]);
    legend(legend_labels, 'Interpreter', 'latex', 'Location','northeast');
    xlabel('$\theta$', 'Interpreter', 'latex');
    ylabel('$\log F(\theta)$', 'Interpreter', 'latex');
    title(['$\epsilon$' sprintf( '= %.2f', eps_DP_vec(i))], 'Interpreter', 'Latex');
    filenametoprint = sprintf('Uniform_width_A_%d_n_%d_eps_%02d', A, n, 100*eps_DP_vec(i));
    print(gcf,'-depsc2', filenametoprint);

    
end

legend(legend_labels, 'Interpreter', 'latex');

%filename = [sprintf('n_%d_param_%s_mu_x_%02d_var_x_%02d_A_%d_FIM_CLT_unif', n, experiment, 10*mu_x, 10*var_x, A) date];
%save(filename);