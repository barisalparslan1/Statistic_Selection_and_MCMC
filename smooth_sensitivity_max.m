function S = smooth_sensitivity_max(X, beta, A_min, A_max)

[N, n] = size(X);

Y = A_max - X;
B_max = A_max - A_min;

Y_sort = [sort(Y, 2) ones(N, 1)*B_max];

temp_1 = log(Y_sort) - (0:n)*beta;
temp_2 = log(Y_sort(:, [2:n+1 n+1]) - Y_sort(:, 1)) - (0:n)*beta;

S = exp(max([temp_1 temp_2], [], 2));