function [theta_vec] = sample_test_theta(y,K,theta,var_theta,prior_mean, ...
    prior_var, a, A, eps_DP, n)

theta_vec = zeros(1,K);

for i = 1:K
    
    theta_prop = theta + sqrt(var_theta)*randn;
    if theta_prop > 0
        if prior_var == inf
            log_prior = 0;
            log_prior_prop = 0;
        else
            log_prior = -0.5*log(2*pi*prior_var) - 0.5*(theta - prior_mean)^2/prior_var;
            log_prior_prop =  -0.5*log(2*pi*prior_var) - 0.5*(theta_prop - prior_mean)^2/prior_var;
        end
        
        mu = theta^a/(a+1);
        mu_prop = theta_prop^a/(a+1);
        
        U = theta^(2*a)*a^2/((a+1)^2*(2*a+1)); 
        U_prop = theta_prop^(2*a)*a^2/((a+1)^2*(2*a+1));
        
        H = U/n + A^(2*a)/(n^2*eps_DP^2);
        H_prop = U_prop/n + A^(2*a)/(n^2*eps_DP^2);
    
        log_likelihood = - 0.5*log(2*pi*H) - 0.5.*(y - mu)^2/H;
        log_likelihood_prop = - 0.5*log(2*pi*H_prop) - 0.5*(y - mu_prop)^2/H_prop;
           
        log_r = (log_likelihood_prop + log_prior_prop) - (log_likelihood + log_prior);
        
        decision = rand < exp(log_r);
        
        if decision == 1
            theta = theta_prop;
        end
    end
    
    theta_vec(i) = theta;
    
end