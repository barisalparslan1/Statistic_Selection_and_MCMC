% FIM approximation and MSE calculation for Section 5.5 in the manuscript
% 
% The population distribution is normal. Depending on the choice, either mean or
% variance can be unknown.
%
% This experiment regards Laplace mechanism with sequntial release of individual  
% data points.
%
% We compare |x|^p for different p as statistics in terms of FIM and MSE.
% 
% The MCMC algorithm is MHAAR-RB.
%
% Use the first code block for obtaining MSE values in Section 5.5.
% Use the second code block for obtaining FIM values in Section 5.5.
% USe the third code block for plotting FIM and MSE together.
%
% Last update: 22 March 2022

clc; clear; close all; fc = 0;
%% MSE Calculation for variance
n = 100;
experiment = 'var'; 
mu_x = 0; var_x  = 2; A = 10;
DP_eps_vec = (0.5:0.5:5); p_vec = [1 2]; 
L_e = length(DP_eps_vec); L_p = length(p_vec);

M = 10000; t_burn = M/4;
N = 100;
theta0 = 10;
sigma_q = 3/sqrt(n);

Num_of_MC_runs = 100;

legends = cell(1, L_p);
Theta_est= zeros(L_p, L_e, Num_of_MC_runs);
E = zeros(L_p, L_e, Num_of_MC_runs);
MSE = zeros(L_e,L_p);

for i = 1:L_p
    p = p_vec(i);
    for j = 1:L_e
        DP_eps = DP_eps_vec(j);

        % compute the sensitivity
        if strcmp(experiment, 'mean') == 1
            Delta = 2*A^p; lap_par_DP = Delta/DP_eps;
        elseif strcmp(experiment, 'var') == 1
            Delta = A^p; lap_par_DP = Delta/DP_eps;
        end

        % run the experiments
        for k = 1:Num_of_MC_runs
            disp([i, j, k]);
            if strcmp(experiment, 'mean') == 1
                X = sqrt(var_x)*randn(1, n) + mu_x;
                y = sign(X).*abs(X).^p + laprnd(1, n, 0, lap_par_DP);
                [Thetas] = MHAAR_RB_DP_norm_mean(y, theta0, DP_eps, A, p, sigma_q, M, N);
                Theta_est(i, j, k) = mean(Thetas(t_burn+1:M));
                E(i, j, k) = Theta_est(i, j, k) - mu_x;
            elseif strcmp(experiment, 'var') == 1
                X = sqrt(var_x)*randn(1, n) + mu_x;
                y = abs(X).^p + laprnd(1, n, 0, lap_par_DP);
                [Thetas] = MHAAR_RB_DP_norm_var(y, theta0, DP_eps, A, p, sigma_q, M, N);
                Theta_est(i, j, k) = mean(Thetas(t_burn+1:M));
                E(i, j, k) = Theta_est(i, j, k) - var_x;
            end
        end
    end
end
MSE(:,1) = mean(squeeze(E(1,:,:)).^2,2);
MSE(:,2) = mean(squeeze(E(2,:,:)).^2,2);

filename = [sprintf('n_%d_param_%s_mu_x_%02d_var_x_%02d_A_%d_MHAAR_RB', n, experiment, 10*mu_x, 10*var_x, A) 'Sequential_MSE' '_' date date];
save(filename);

%% FIM Calculation 
fc = 0;
theta = 2; %true theta 

A = 10;
n = 100;  p_vec = [1 2]; L_p = length(p_vec);
DP_eps_vec = (0.5:0.5:5); L_e = length(DP_eps_vec);

% number of proposals for the latent variable
N = 1000;
% number of MC runs
M = 1000;

FIM_V = zeros(L_p,L_e);

for i = 1:L_p
    p = p_vec(i);
    for j = 1:L_e
        disp([i, j]);
        DP_eps = DP_eps_vec(j);
        
        %Calculating FIM values using MCMC algorithm
        FIM = FIM_DP_Lap_seq_release(theta, p, n, A, DP_eps, N, M);
        FIM_V(i,j) = FIM ;
    end
end

filename = [sprintf('theta_%02d_A_%d',theta, A) '_sequential_FIM' '_' date];
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
            plot(log(MSE(:, i)), '-');
            hold on;
        end
        hold off;
        set(gca, 'Xtick', 1:2:L_e, 'XtickLabel', DP_eps_vec(1:2:L_e));
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('$\log$(MSE)', 'Interpreter','Latex')
        legend(legends,'Interpreter', 'latex');
    else
        for i = 1:L_p
            plot(log(FIM_V(i, :)), '-');
            hold on;
        end
        hold off
        set(gca, 'Xtick', 1:2:L_e, 'XtickLabel', DP_eps_vec(1:2:L_e));
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('$\log(F(\theta))$', 'Interpreter','Latex')
        legend(legends,'Interpreter', 'latex','Location','northwest'); 
    end
end
