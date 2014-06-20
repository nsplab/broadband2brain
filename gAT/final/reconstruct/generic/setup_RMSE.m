function [setup_obj] = setup_RMSE(T, window, args)
% Precomputes values for generic kernel
% args.h is function handler for kernel
% args.h_args is args for kernel

% copy/pasted from setup_ker.m
% not all parts of this are actually used - this is kind of a hack

L = args.L;
h = args.h;
t_resolution = args.t_resolution;

% Set up for convolution
convolve_setup(T, window, L, h, args);

T_s = T / L;  % Spacing of samples

% Covariance matrix for samples y
% Samples assumed to follow a multivariate gaussian distribution
cov_inv = eye(L);  % changed (for RMSE)

% Compute mu, used to correct for non-zero-mean noise in samples
mu = ones(L, 1) / sqrt(2 * pi);

% Discrete values of t
t_sampled = 0 : t_resolution : T;
num_t_bins = length(t_sampled);

% Precompute values for grid search to estimate t_k
B_full = zeros(L, num_t_bins);
for n = 1 : L,
    for i = 1 : num_t_bins,
        B_full(n, i) = h(n*T_s - t_sampled(i), args);
    end
end
B_const = zeros(1, num_t_bins);
for i = 1 : num_t_bins,
    B_const(1, i) = B_full(:, i)' * cov_inv * B_full(:, i);
end

% Precompute values of kernel
h_mat = zeros(length(t_sampled), L);
for i = 1 : length(t_sampled)
    for n = 1 : L
        h_mat(i, n) = h(n*T_s - t_sampled(i), args);
    end
end

setup_obj = struct();
setup_obj.cov_inv = cov_inv;
setup_obj.mu = mu;
setup_obj.t_sampled = t_sampled;
setup_obj.B_full = B_full;
setup_obj.B_const = B_const;
setup_obj.h_mat = h_mat;

end