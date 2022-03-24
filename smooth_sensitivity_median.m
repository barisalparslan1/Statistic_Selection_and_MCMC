function S = smooth_sensitivity_median(X, beta, A_min, A_max)

[N, n] = size(X);
X = sort(X, 2);

X = X - A_min;

X = [zeros(N, 1) X ones(N, 1)*(A_max-A_min)];
m = floor((n+3)/2);

temp = zeros(N, n+1);

for k = 0:n
    t_vec = 0:(k+1);        
    temp(:, k+1) = max(X(:, min(n+2, m+t_vec)) - X(:, max(1, m+t_vec-k-1)), [], 2);
end

S = exp(max(log(temp)-beta*(0:n), [], 2));