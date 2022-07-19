% strategy comparison / new experiment

theta_vals = 0.1:0.1:1; L_e = length(theta_vals);
A = 10;
eps_DP = 1;
a_vec = 0.1:0.1:2; L_a = length(a_vec);
n = 100;
n1 = n*0.1;
n2 = n*0.9;
M = 15000;

% for MC run
theta0 = 2;
var_theta = 10/sqrt(n);
K = 10000;
t_burn = K/4;

theta_est_s1 = zeros(1,M);
theta_est_s2 = zeros(1,M);
t_burn_post = K/4;

Errors_s1 = zeros(1,M);
Errors_s2 = zeros(1,M);
MSE = zeros(2,L_e);

c1 = 0; c2 = 0;

temp_count_vec = zeros(M,L_e);
for k = 1: L_e
    theta = theta_vals(k);
    disp(theta);
    % strategy 1
    for i = 1:M
        X = -theta + (theta+theta)*rand(1,n);
        X1 = X(1:n1);
        X2 = X(n1+1:end);

        Delta = A/n1; 
        v = randn*sqrt(Delta^2/eps_DP^2); 
        y1 = mean(abs(X1)) + v;

        [temp_vec] = sample_test_theta(y1,K,theta0,var_theta, 0, inf, 1, A,eps_DP,n1);
        theta_vec = temp_vec(t_burn+1:end);
        L_theta = length(theta_vec);

        sample_mean = mean(theta_vec);
        sample_var = var(theta_vec);
        
        temp_sum_vec = zeros(1,L_a);
        for j = 1:L_theta
            theta_temp = theta_vec(j);
            mu_deriv = a_vec.*(theta_temp).^(a_vec-1)./(a_vec+1);
            U = theta_temp.^(2*a_vec).*(a_vec.^2./((a_vec+1).^2.*(2*a_vec + 1)))/n2 + A.^(2*a_vec)/(n2^2*eps_DP^2);
            U_deriv = theta_temp.^(2*a_vec-1).*(2*a_vec.^3./((a_vec+1).^2.*(2*a_vec + 1)))/n2;
            FIM_temp = mu_deriv.^2./U + 0.5*(U_deriv.^2)./(U.^2);

            for p=1:length(FIM_temp)
                temp_sum_vec(p) = temp_sum_vec(p)+1/FIM_temp(p);
            end
        end
        
        [argvalue, argmin] = min(temp_sum_vec);
        
        a_val = a_vec(argmin);
        temp_count_vec(i,k) = a_val;

        Delta = A^a_val/n2; 
        v = randn*sqrt(Delta^2/eps_DP^2);
        y2 = mean(abs(X2).^a_val) + v;
        [theta_post] = sample_test_theta(y2,K,theta0,var_theta,sample_mean,sample_var,a_val,A,eps_DP,n2);

        theta_post_vec = theta_post(t_burn_post+1:end);
        theta_est_s1(i) = mean(theta_post_vec);

        Errors_s1(i) = theta_est_s1(i) - theta; 
    end

    % strategy 2 
    for i = 1:M
        X = -theta + (theta+theta)*rand(1,n);
        Delta = A/n;
        v = randn*sqrt(Delta^2/eps_DP^2);
        y = mean(abs(X)) + v;

        [temp_vec] = sample_test_theta(y,K,theta0,var_theta,0, inf, 1, A,eps_DP,n);

        theta_vec = temp_vec(t_burn+1:end);

        theta_post_vec = theta_vec(t_burn_post+1:end);
        theta_est_s2(i) = mean(theta_post_vec);

        Errors_s2(i) = theta_est_s2(i) - theta; 
    end

    MSE(1,k) = mean(Errors_s1.^2);
    MSE(2,k) = mean(Errors_s2.^2);
end

figure(1);
subplot(1, 2, 1);
for i = 1:2
   plot(theta_vals,log(MSE(i, :)), '*-');
   hold on;
end
hold off  
xlabel('true value of $\theta$', 'Interpreter','Latex');
ylabel('$\log$(MSE)', 'Interpreter','Latex')
legend('stat selection','no stat selection','Location','southeast');

subplot(1, 2, 2);
boxplot(temp_count_vec,theta_vals, 'symbol', '')
xlabel('true value of $\theta$', 'Interpreter','Latex');
ylabel('selected $a$', 'Interpreter','Latex')
set(gca, 'ylim',[0.25, 1.25])
filename = ('theta0.1_1_a0.1_2_M_15000_K_10000_strategy_comparison.mat');
save(filename);
