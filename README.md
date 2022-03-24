# Statistic Selection and MCMC for Differentially Private Bayesian Estimation

# Brief information about the work
In this work we focused on Bayesian estimation of an unknown population parameter when a certain statistic is shared with a noise to provide differential privacy. We compared statistics in terms of Fisher information and Mean Squared Error of their Bayesian estimator. For estimating numerical values of Fisher information and  parameters, we utilized state-of-art MCMC algorithms. For details and the results please check the manuscript.  

Authors: Barış Alparslan and Sinan Yıldırım (Faculty of Engineering and Natural Sciences, Sabancı University, Istanbul, Turkey).

# main files

#### main_CLT_norm_mean.m for the calculations in Example 1:
FIM calculation for Example 1 in the manuscript. The MCMC algorithm is a MH algorithm.

- Use the first code block for plots in Example 1 in the manuscript.

#### main_CLT_norm_var for the experiments in Section 5.1:
FIM and MSE calculation for Example 2 and Section 5.1 in the manuscript. The MCMC algorithm is a MH algorithm.

- Use the first code block for obtaining FIM values for Example 2
- Use the second code block for obtaining MSE values in Section 5.1
- Use the third code block for plotting FIM and MSE values together

#### main_FIM_CLT_unif.m for the experiments in Example 3:

- Use the first code block for obtaining FIM values for Example 3

#### main_FIM_binary.m for the experiments in Example 4:
This is the experiment regarding binary responses.

- Use the first code block for obtaining FIM values for Example 4.

#### main_CLT_Lap_norm_var.m for the experiments in Section 5.2 and Section 5.3:
FIM approximation, MSE calculation and IAC calculation for Section 5.2 and 5.3 in the manuscript. The MCMC algorithms are PMMH and MHAAR.

- Use the first code block for obtaining FIM values in Section 5.2
- Use the second code block for obtaining MSE values in Section 5.2
- Use the third code block for IAC values at the Section 5.3.
- Use the code block at the bottom for plotting FIM and MSE together.

#### main_non_additive_norm_var.m for the experiments in Section 5.4:
FIM approximation, MSE and ACF calculation for Section 5.4 in the manuscript. The MCMC algorithm is MHAAR.

- Use the first code block for obtaining FIM values in Section 5.4
- Use the second code block for obtaining MSE values in Section 5.4
- Use the third code block for ACF values at the Section 5.4 with median.
- Use the code block at the bottom for plotting FIM and ACF together.

#### main_MHAAR_RB_DP.m for the experiments in Section 5.5:
FIM approximation and MSE calculation for Section 5.5 in the manuscript. The MCMC algorithm is MHAAR-RB.

- Use the first code block for obtaining MSE values in Section 5.5.
- Use the second code block for obtaining FIM values in Section 5.5.
- USe the third code block for plotting FIM and MSE together.

# functions 
## functions for FIM calculation
#### FIM_DP_norm_mean_CLT 
###### Inputs: DP level, data size, moment of x, boundary of x, mean parameter
###### output: Fisher information

This function calculates the FIM of the noisy summary statistic(additive) with respect to the mean parameter of the normal distribution with variance 1. Noise: Gaussian

#### FIM_DP_norm_var_CLT 
###### Inputs: DP level, data size, moment of x, boundary of x, variance parameter
###### output: Fisher information

This function calculates the FIM of the noisy summary statistic(additive) with respect to the variance parameter of the normal distribution with mean 0. Noise: Gaussian

#### FIM_DP_CLT_Lap_norm_var
###### Inputs: true variance parameter, moment of x, data size, boundary of x, DP level, size of latent variable, number of MCMC runs
###### output: Fisher information

This function calculates the FIM of the noisy summary statistic(additive) with respect to the variance parameter of the normal distribution with mean 0. Noise: Laplace

#### FIM_DP_norm_var
###### Inputs: true variance parameter, moment of x, data size, boundary of x, DP level, DP variable for smooth sensitivity, size of latent variable, number of MCMC runs, name of the statistic
###### output: Fisher information

This function calculates the FIM of the noisy individual statistic(non-additive) with respect to the variance parameter of the normal distribution with mean 0. Noise: Laplace

#### FIM_DP_Lap_seq_release
###### Inputs: true variance parameter, moment of x, data size, boundary of x, DP level, size of latent variable, number of MCMC runs
###### output: Fisher information

This function calculates the FIM of the noisy individual statistic(sequential release) with respect to the variance parameter of the normal distribution with mean 0. Noise: Laplace

## functions for Bayesian estimation and MSE calculation
#### MH_DP_CLT_norm_var
###### Inputs: shared statistic, initial theta, moment of x, data size, boundary of x, DP level, number of MCMC runs, variance for proposal distribution
###### output: Fisher information

This function calculates the MSE of the estimation for the variance parameter of the normal distribution when the summary statistic(additive) is shared with Gaussian noise. MCMC: Metropolis-Hastings

#### PMMH_DP_CLT_Lap_norm_var
###### Inputs: shared statistic, initial theta, moment of x, data size, boundary of x, DP level, number of MCMC runs, variance for proposal distribution, number of latent variables
###### output: Fisher information

This function calculates the MSE of the estimation for the variance parameter of the normal distribution when the summary statistic(additive) is shared with Laplace noise. MCMC: Pseudo-marginal Metropolis-Hastings

#### MHAAR_DP_CLT_Lap_norm_var
###### Inputs: shared statistic, initial theta, moment of x, data size, boundary of x, DP level, number of MCMC runs, variance for proposal distribution, number of latent variables
###### output: Fisher information

This function calculates the MSE of the estimation for the variance parameter of the normal distribution when the summary statistic(additive) is shared with Laplace noise. MCMC: Metropolis-Hastings with Averaged Acceptance Ratio

#### DM_MHAAR_norm_var
###### Inputs: shared statistic, moment of x, boundary of x, data size, number of MCMC runs, number of latent variables, DP level, DP variable for smooth sensitivity, initial theta,variance for proposal distribution, name of the statistic
###### output: Fisher information

This function calculates the MSE of the estimation for the variance parameter of the normal distribution when the summary statistic(non-additive) is shared with Laplace noise. MCMC: Metropolis-Hastings with Averaged Acceptance Ratio
