% FIM approximation, MSE calculation and IAC calculation for Section 5.2
% and 5.3 in the manuscript
% 
% The population distribution is normal with zero mean and unknown variance theta.
% Laplace mechanism is used to provide DP.
% 
% This experiment regards Laplace mechanism with additive statistics. 
% 
% We compare |x|^p for different p as statistics in terms of FIM and MSE.
% 
% The MCMC algorithm are PMMH and MHAAR.
%
% Use the first code block for obtaining FIM values in Section 5.2
% Use the second code block for obtaining MSE values in Section 5.2
% Use the third code block for IAC values at the Section 5.3.
% Use the code block at the bottom for plotting FIM and MSE together.
%
% Last update: 22 March 2022

clear; clc; close all; fc = 0;

%%  Experiments for FIM calculations

theta = 2; %true theta 

A = 10;
n = 100;  p_vec = [1 2]; L_p = length(p_vec);
DP_eps_vec = [0.25:0.25:5, 5.5:0.5:10]; L_e = length(DP_eps_vec);

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
        FIM = FIM_DP_CLT_Lap_norm_var(theta, p, n, A, DP_eps, N, M);
        FIM_V(i,j) = FIM ;
    end
end

filename = [sprintf('n_%d_var_x_%02d_A_%d_eps_%d', n, theta, A,DP_eps) 'Lap_additive_FIM' '_' date];
save(filename);
%%  Experiments for MSE calculations for MCMC estimates
clear; clc; close all; fc = 0;

theta = 2; %true theta 

A = 10;
n = 100;  p_vec = [1 2]; L_p = length(p_vec);
DP_eps_vec = [0.25:0.25:5, 5.5:0.5:10]; L_e = length(DP_eps_vec);

%Sigma for proposal theta generation
sigma_q = 11/sqrt(n);

% Number of MCMC iterations
K = 10000; t_burn = K/4;

% number of proposals for the latent variable
N = 100;
% number of MC runs (to average)
M = 100;

%initializations
theta0 = 2;
mu0 = 1;  %mu0 is input for MHAAR_DP


Acc = zeros(L_p, L_e, M);
MSE = zeros(L_p, L_e);
Errors = zeros(L_p, L_e, M);

for i = 1:L_p
    p = p_vec(i);
    for j = 1:L_e
        disp([i, j]);
        DP_eps = DP_eps_vec(j);
        for m = 1:M
            %generating noisy observations
            X = sqrt(theta)*randn(1,n);
            Delta = A^p/n;
            v = laprnd(1,1,0,Delta/DP_eps);
            y = mean(abs(X).^p) + v;
            
            %Run MCMC algorithm - PMMH_DP for calculating MSEs
            [outputs] = PMMH_DP_CLT_Lap_norm_var(y, theta0, p, n, A, DP_eps,K,sigma_q,N);
            
            %Run MCMC algorithm - MHAAR_DP for calculating MSEs
            %[outputs] = MHAAR_DP_CLT_Lap_norm_var(y, theta0,mu0, p, n, A, DP_eps(dp),K,sigma_q,N);
            
            %Estimate parameter and corresponding error
            Thetas_conv = outputs.Thetas(t_burn + 1:end);
            theta_est = mean(Thetas_conv);
            Errors(i, j, m) = theta_est - theta;
            Acc(i,j,m) = mean(outputs.dec_vec(t_burn + 1:end));
        end
        MSE(i, j) = mean(Errors(i, j,:).^2);
    end
end

filename = [sprintf('n_%d_var_x_%02d_A_%d_eps_%d', n, theta, A,DP_eps) 'Lap_additive_MSE' '_' date];
save(filename);
%% IAC Score comparison for PMMH_DP and  MHAAR_DP for normal variance

clear; clc; close all; fc = 0;

theta = 2; %true theta 

A = 10;
n = 100;  p = 1; 
DP_eps = 5;
Delta = A^p/n;

%Sigma for proposal theta generation
sigma_q = 20/sqrt(n);

%Number of MCMC iterations
K = 200000; t_burn = K/4;

% number of MC runs (to average)
M = 100;
%Number of proposals for latent variable
N = [2, 5, 10, 20, 50, 100]; L_N = length(N);

%initializations
theta0 = 2;
mu0 = 1;

tau_s = zeros(L_N,2);
for t=1:L_N
    disp(N(t))
    
    %generating noisy observations
    X = sqrt(theta)*randn(1,n);
    v = laprnd(1,1,0,Delta/DP_eps);
    y = mean(abs(X).^p) + v;
    
    %Run MCMC algorithms with corresponding N and y 
    [outputs_A2] = PMMH_DP_CLT_Lap_norm_var(y, theta0, p, n, A, DP_eps,K,sigma_q,N(t)); 
    [outputs_A3] = MHAAR_DP_CLT_Lap_norm_var(y, theta0,mu0, p, n, A, DP_eps,K,sigma_q,N(t)); 
    
    %Estimate the parameter
    Thetas_conv_A3 = outputs_A3.Thetas(t_burn + 1:end); Thetas_conv_A2 = outputs_A2.Thetas(t_burn + 1:end);
    
    %Calculate IAC scores for two MCMC algorithms 
    [tau_A2, ~, cf_A2, ~] = IAC_Sokal(transpose(Thetas_conv_A2), 6, 'unbiased');
    [tau_A3, ~, cf_A3, ~] = IAC_Sokal(transpose(Thetas_conv_A3), 6, 'unbiased');
    tau_s(t,1)=tau_A2; tau_s(t,2)=tau_A3;
end

filename = [sprintf('n_%d_var_x_%02d_A_%d_eps_%d', n, theta, A,DP_eps) 'Lap_additive_IAC' '_' date];
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
            plot(log(MSE(i, :)), '-');
            hold on;
        end
        hold off  
        set(gca, 'Xtick', [2:6:L_e,L_e], 'XtickLabel', DP_eps_vec([2:6:L_e,L_e]));
        set(gca,'ylim',[-3,2.5]);
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('$\log$(MSE)', 'Interpreter','Latex')
        legend(legends,'Interpreter', 'latex');
    else
        for i = 1:L_p
            plot(log(FIM_V(i, :)), '-');
            hold on;
        end
        hold off
        set(gca, 'Xtick', [2:6:L_e,L_e], 'XtickLabel', DP_eps_vec([2:6:L_e,L_e]));
        set(gca,'ylim',[-2,2.5]);
        xlabel('$\epsilon$', 'Interpreter','Latex');
        ylabel('$\log(F(\theta))$', 'Interpreter','Latex')
        legend(legends,'Interpreter', 'latex','Location','southeast'); 
    end
end