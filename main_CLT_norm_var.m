% FIM and MSE calculation for Example 2 and Section 5.1 in the manuscript 
%
% The population distribution is normal with zero mean and unknonw variance theta. 
% Gaussian mechanism is
% used to provide eps, delta -DP.
%
% We compare |x|^p for different p as statistics in terms of FIM and MSE.
%
% The MCMC algorithm is a MH algorithm.
%
% Use the first code block for obtaining FIM values for Example 2
% Use the second code block for obtaining MSE values in Section 5.1
% Use the third code block for plotting FIM and MSE values together
%
% Last update: 22 March 2022

clear; clc; close all; fc = 0;
%% FIM values in Example 2 and Section 5.1
A = 10;
n = 100; p_vec = [1 2 3]; L_p = length(p_vec);
eps_DP_vec = [1 inf]; L_eps = length(eps_DP_vec);
legends = cell(1, L_p);

Theta = 0.1:0.1:10; L_theta = length(Theta);

FIM = zeros(L_eps, L_p, L_theta);

for i1 = 1:L_eps
    eps_DP = eps_DP_vec(i1);
    for i2 = 1:L_p
        p = p_vec(i2);
        for i3 = 1:L_theta
            theta = Theta(i3);
            FIM(i1, i2, i3) = FIM_DP_norm_var_CLT(theta, n, eps_DP, p, A);
        end
    end
end

% To obtain FIM values in Section 5.1 change Theta as 2 and extend
% eps_DP_vec
% 
%
filename = [sprintf('n_%d_A_%d', n,A) 'CLT_norm_var_FIM' '_' date];
save(filename);

%% MSE values in Section 5.1
clear; clc; close all; fc = 0; 

theta = 2; %true theta 

A = 10;
n = 100; p_vec = [1 2]; L_p = length(p_vec);
DP_eps_vec = [0.25:0.25:5, 5.5:0.5:10 inf]; L_e = length(DP_eps_vec);

%Sigma for proposal theta generation
sigma_q = 10/sqrt(n);

% Number of MCMC iterations
K = 100000; t_burn = K/4;
% number of MC runs
M = 1000;

%initialization
theta0 = 2;
MSE = zeros(L_p, L_e);

legends = cell(1, L_p);
Errors = zeros(L_p, L_e, M);
Acc = zeros(L_p, L_e, M);

for i = 1:L_p
    p = p_vec(i);
    legends{i} = sprintf('p = %d', p);
    for j = 1:L_e
        disp([i, j])
        DP_eps = DP_eps_vec(j);
        for m = 1:M
            %noisy observation generation
            X = sqrt(theta)*randn(1,n);
            Delta = A^p/n;
            v = randn*sqrt(Delta^2/DP_eps^2);
            y = mean(abs(X).^p) + v;
            
            %Run MCMC algorithm
            [outputs] = MH_DP_CLT_norm_var(y, theta0, p, n, A, DP_eps, K, sigma_q);
            
            %Estimate parameter and MSE
            Thetas_conv = outputs.Thetas(t_burn + 1:end);
            theta_est = mean(Thetas_conv);
            Errors(i, j, m) = theta_est - theta;
            Acc(i, j, m) = mean(outputs.dec_vec(t_burn + 1:end));
        end
        MSE(i, j) = mean(Errors(i, j, :).^2);
    end
end

filename = [sprintf('n_%d_A_%d_theta_%d_eps_%d_', n, A, theta, DP_eps) 'CLT_norm_var_MSE' '_' date];
save(filename);

%% Plot FIM and MSE together

for i = 1:L_p
    p = p_vec(i); 
    legends{i} = sprintf('$a$ = %d', p);
end
for k = 1:2
    subplot(1,2,k);
    if k == 1
        for i = 1:L_p
            plot((MSE(i, :)), '-');
            hold on;
        end
        hold off  
        set(gca, 'Xtick', [2:6:L_e,L_e], 'XtickLabel', DP_eps_vec([2:6:L_e,L_e]));
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('MSE')
        legend(legends,'Interpreter', 'latex');
    else
        for i = 1:L_p
            plot(log(FIM(:,i,1)), '-');
            hold on;
        end
        hold off
        set(gca, 'Xtick', [2:6:L_e,L_e], 'XtickLabel', eps_DP_vec([2:6:L_e,L_e]));
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('$\log(F(\theta))$', 'Interpreter','Latex')
        legend(legends,'Interpreter', 'latex','Location','southeast'); 
    end
end
