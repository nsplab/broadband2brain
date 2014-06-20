function [cov_inv, mu, t_sampled, B_full, B_const] = FRI_sinc_diode_setup(L, T, B, t_resolution)
% Precomputes values for use by FRI_sinc_diode.m

% Inputs:
% (most same as in FRI_sinc_diode.m)
% t_resolution - separation between possible discrete t_k values

T_s = T / L;  % Spacing of samples

% Covariance matrix for samples y
% Samples assumed to follow a multivariate gaussian distribution
cov = zeros(L, L);
for i = 1 : L,
    for j = 1 : L,
        cov(i, j) = normalized_sinc(abs(i-j)*T_s, B);
    end
end

cov_inv = cov^(-1);
%cov_inv = eye(size(cov_inv));  % TESTING

% Compute mu, used to correct for non-zero-mean noise in samples
mu = ones(L, 1) / sqrt(2 * pi);

% Discrete values of t
t_sampled = 0 : t_resolution : T;
num_t_bins = length(t_sampled);

% Precompute values for grid search to estimate t_k
B_full = zeros(L, num_t_bins);
for n = 1 : L,
    for i = 1 : num_t_bins,
        B_full(n, i) = normalized_sinc((n-1/2)*T_s - t_sampled(i), B);
    end
end
B_const = zeros(1, num_t_bins);
for i = 1 : num_t_bins,
    B_const(1, i) = B_full(:, i)' * cov_inv * B_full(:, i);
end

end
