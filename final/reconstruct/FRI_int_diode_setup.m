function [cov_inv, mu, t_sampled, B_full, B_const] = FRI_sinc_diode_setup(L, T, t_resolution)
% Precomputes values for use by FRI_int_diode.m

% Covariance matrix for samples y
% Samples assumed to follow a multivariate gaussian distribution
cov = zeros(L, L);
for i = 0 : L-1,
    for j = 0 : L-1,
        cov(i+1, j+1) = T^(i+j) / (factorial(i) * factorial(j) * (i+j+1));
    end
end

cov_inv = cov^(-1);

% Compute mu, used to correct for non-zero-mean noise in samples
mu = zeros(L, 1);
for l = 0 : L-1,
    mu(l+1) = T^(l+1)/(sqrt(2*pi)*factorial(l+1));
end

% Discrete values of t
t_sampled = 0 : t_resolution : T;
num_t_bins = length(t_sampled);

% Precompute values for grid search to estimate t_k
B_full = zeros(L, num_t_bins);
for l = 0 : L-1,
    for i = 1 : num_t_bins,
        B_full(l+1, i) = 1/factorial(l) * (T - t_sampled(i))^l;
    end
end
B_const = zeros(1, num_t_bins);
for i = 1 : num_t_bins,
    B_const(1, i) = B_full(:, i)' * cov_inv * B_full(:, i);
end

end
