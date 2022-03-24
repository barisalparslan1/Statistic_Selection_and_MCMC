% FIM approximation, MSE and ACF calculation for Section 5.4 in the manuscript
% 
% The population is normal with zero mean and unknown variance theta.
% Laplace mechanism is used to provide DP.
% 
% This experiment regards Laplace mechanism with non-additive statistics using  
% smooth sensitivity
%
% We compare |x|^p for different p as statistics in terms of FIM and MSE.
% 
% The MCMC algorithm is MHAAR.
%
% Use the first code block for obtaining FIM values in Section 5.4
% Use the second code block for obtaining MSE values in Section 5.4
% Use the third code block for ACF values at the Section 5.4 with median.
% Use the code block at the bottom for plotting FIM and ACF together.
%
% Last update: 22 March 2022

clear; clc; close all; fc = 0;
%% Experiments for FIM calculations
A = 10;
n = 100;
DP_eps_vec = 5; L_e = length(DP_eps_vec);
DP_delta = 1/n^2;

Theta_vec = 0.25:0.25:5; L_t = length(Theta_vec);
p_vec = 1; L_p = length(p_vec);

legends = {'median' 'max'};

% number of proposals for the latent variable
N = 200;

% number of MC runs (to average)
M = 1000;

FIM_V = zeros(2, L_t);

for i1 = 1:L_p
    p = p_vec(i1);
    for i2 = 1:L_e
        DP_eps = DP_eps_vec(i2);
        for i3 = 1:L_t
            disp([i1, i2, i3]);
            theta = Theta_vec(i3);
            FIM_V(1, i3) = FIM_DP_norm_var(theta, p, n, A, DP_eps, DP_delta, N, M, 'median');
            FIM_V(2, i3) = FIM_DP_norm_var(theta, p, n, A, DP_eps, DP_delta, N, M, 'max');
        end
    end
end

% fc = fc + 1; figure(fc);
% plot(Theta_vec, log(FIM_V(1, :)), '.-', Theta_vec, log(FIM_V(2, :)), '*-');
% 
% xlabel('$\theta$', 'Interpreter','Latex');
% ylabel('$\log(F(\theta))$', 'Interpreter','Latex')
% title('FIM values for different thetas');
% legend(legends);

filename = [sprintf('n_%d_%02d_A_%d_', n, A) 'Lap_non_additive_FIM' '_'  date];
save(filename);
%% Experiments for MSE calculations for MCMC estimates
theta_true = 2;
Stats_names = {'median', 'max'}; L_s = length(Stats_names);

A = 10;
n = 100; p_vec = 1; L_p = length(p_vec);
DP_eps_vec = 5; L_e = length(DP_eps_vec);
DP_delta = 1/(n^2);

%Sigma for proposal theta generation
sigma_q = 20/sqrt(n);

% Number of MCMC iterations
K = 10000; t_burn = K/4;
% number of proposals for the latent variable
N = 50;

% initial value
theta0 = theta_true;

% number of MC runs (to average)
M = 100;

Errors = zeros(L_p, L_e, L_s, M);
Acc = zeros(L_p, L_e, L_s, M);
tau= zeros(L_p, L_e, L_s, M);
cf= cell(L_p, L_e, L_s, M);


for i1 = 1:L_p
    p = p_vec(i1);
    A_sens = A^p;

    for i2 = 1:L_e
        DP_eps = DP_eps_vec(i2);

        beta = DP_eps/(2*log(2/DP_delta));
        alfa = DP_eps/2;

        for i3 = 1:L_s
            stats_name = Stats_names{i3};

            for i4 = 1:M
                disp([i1, i2, i3, i4]);

                X = sqrt(theta_true)*randn(1, n);
                X = min(A, max(-A, X));
                S = abs(X).^p;

                %Generate noisy observations using smooth sensitivity
                if strcmp(stats_name, 'median') == 1
                    sigma = smooth_sensitivity_median(S, beta, 0, A_sens);
                    v = laprnd(1, 1, 0, sigma/alfa);
                    y = median(S) + v;
                elseif strcmp(stats_name, 'max') == 1
                    sigma = smooth_sensitivity_max(S, beta, 0, A_sens);
                    v = laprnd(1, 1, 0, sigma/alfa);
                    y = max(S) + v;
                end

                % run the MCMC algorithm
                [outputs] = DM_MHAAR_norm_var(y, p, A, n, K, N, DP_eps, ...
                    DP_delta, theta0, sigma_q, stats_name);

                %Parameter estimation
                Thetas_conv = outputs.Thetas(t_burn + 1:end);
                theta_est = mean(Thetas_conv);

                Errors(i1, i2, i3, i4) = theta_est - theta_true;
                Acc(i1, i2, i3, i4) = mean(outputs.dec_vec(t_burn + 1:end));
                
                [tau(i1, i2, i3, i4), ~, cf{i1, i2, i3, i4}, ~] = IAC_Sokal(Thetas_conv, 6, 'unbiased');

            end
        end
    end
end

%MSE calculation for median and max
errors_median = squeeze(Errors(:,:,1,:));
MSE_median = mean(errors_median.^2);

errors_max = squeeze(Errors(:,:,2,:));
MSE_max = mean(errors_max.^2);

filename = [sprintf('n_%d_var_x_%02d_A_%d_', n, theta_true, A) 'Lap_non_additive_MSE' '_' date];
save(filename);

%% ACF calculation for DM_MHAAR for median

theta = 2; %true theta

A = 10;
n = 100;  p = 1;
DP_eps = 5; DP_delta = 1/(n^2);
beta = DP_eps/(2*log(2/DP_delta)); alfa = DP_eps/2;
A_sens = A^p;

%Sigma for proposal theta generation
sigma_q = 20/sqrt(n);

%Number of MCMC iterations
K = 10000; t_burn = K/4;
%Number of proposals for latent variable
N = 50;

% number of MC runs (to average)
M = 5;

stat = {'median', 'max'};

%initial values
theta0 = 2;

cf_f = zeros(M, (3*K/4));
cf_final_score = zeros((3*K/4),2);

for st = 1:length(stat)
    for t=1:M
        disp(t)
        X = sqrt(theta)*randn(1,n);
        S = abs(X).^p;
        
        % compute smooth sensitivity and generate noisy observations
        if strcmp(stat{st}, 'median')== 1
            sigma = smooth_sensitivity_median(S, beta, 0, A_sens);
            v = laprnd(1,1,0,sigma/alfa);
            y = median(abs(X).^p) + v;
        elseif strcmp(stat{st}, 'max') == 1
            sigma = smooth_sensitivity_max(S, beta, 0, A_sens);
            v = laprnd(1,1,0,sigma/alfa);
            y = max(abs(X).^p) + v;
        end
        
        % run the MCMC algorithm
        [outputs] = DM_MHAAR_norm_var(y, p, A, n, K, N, DP_eps, ...
            DP_delta, theta0, sigma_q, stat{st});
        Thetas_conv = outputs.Thetas(t_burn + 1 :end);

        [tau, ~, cf, ~] = IAC_Sokal(Thetas_conv, 6, 'unbiased');

        cf_f(t,:) = cf;
    end

    cf_final_score(:,st) = mean(cf_f,1);
end

filename = [sprintf('n_%d_var_x_%02d_A_%d_eps_%d', n, theta, A,DP_eps) 'Lap_non_additive_IAC' '_' date];
save(filename);


%% Plotting FIM and ACF Together
legends = {'median' 'max'};

for k = 1:2
    subplot(1,2,k);
    if k == 1
        plot(Theta_vec, log(FIM_V(1, :)), '.-', Theta_vec, log(FIM_V(2, :)), '*-');
        set(gca,'ylim',[-3,7])
        xlabel('$\theta$', 'Interpreter','Latex');
        ylabel('$\log(F(\theta))$', 'Interpreter','Latex')
        legend(legends);
    else
        for st = 1:2
            set(gca,'XLim',[0 4000])
            plot(cf_final_score(:,st),'LineWidth',1);
            hold on
        end
        hold off
        set(gca,'ylim',[-0.1,0.5])
        ylabel('ACF');
        legend(legends);
    end
end