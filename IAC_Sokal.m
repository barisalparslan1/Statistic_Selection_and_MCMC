function [tau, v, cf, M] = IAC_Sokal(x, c, norm_type)

% [tau, v, cf, M] = IAC_Sokal(x, c, norm_type)
% 
% This function estimatates the integrated autocorrelation (IAC) time of 
% a WSS process X_1, X_2, ... from a realisation of it, x. The input c 
% determines window size M in the Sokal sense. norm determines the 
% whether unbiased or biased normalisation will be performed in order to 
% estimate the autocorrelation
% 
% tau: IAC time
% v: the estimated variance of X,
% cf: autocorrelation function
% M: window size
% 
% Sinan Yildirim, 13.05.2014
% Last update: 21.05.2014, 21.39

% Find the autocorrelation function of x via fft and ifft operations:

% First get rid of the mean
x = x - mean(x);

% fft with zero padding:
n = length(x);
f = fft(x, 2^nextpow2(2*n-1));

% Find the periodogram
pf = real(f).^2+imag(f).^2;

% find the inverse transform to get the autocorrelation function:
cf = ifft(pf);
% get rid of the zeros padding
cf = cf(1:n);

cf = cf(:);
v = cf(1)/n;


% normalise
if nargin == 2
    norm_type = 'unbiased';
end

if strcmp(norm_type, 'biased') == 1
    % biased estimate of the autocorrelation function
    cf = cf./cf(1);
elseif strcmp(norm_type, 'unbiased') == 1
    % unbiased estimate of the autocorrelation function
    cf = cf./(cf(1)*(n - (0:n-1)')/n);
end

% Find the window size
tau = 1;
m = 0;
check = tau <= 2*(m+1)/c;

while check == 0 && m+1 <= n-1
    m = m + 1;
    tau = tau + 2*cf(m+1);
    check = tau < 2*(m+1)/c;
end
M = m;

