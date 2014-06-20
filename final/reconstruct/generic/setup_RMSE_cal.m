function [setup_obj] = setup_RMSE_cal(T, window, args)
% Precomputes values for generic kernel
% args.h is function handler for kernel
% args.h_args is args for kernel

% copy/pasted from setup_ker.m
% not all parts of this are actually used - this is kind of a hack

L = args.L;
h = args.h;
t_resolution = args.t_resolution;

% Set up for convolution (comment because unneeded and SLOW) ***
% comment this line when using calcium data?
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
%B_full = zeros(L, num_t_bins);
B_full = h(repmat((1 : L)', [1 num_t_bins])*T_s - repmat(t_sampled, [L 1]), args);
%for n = 1 : L,
%    for i = 1 : num_t_bins,
%        B_full(n, i) = h(n*T_s - t_sampled(i), args);
%    end
%end

B_const = zeros(1, num_t_bins);
% TODO: Vectorize
%B_const2 = B_full' .^ 2 * cov_inv;
%size(B_full)
%size(cov_inv)
%num_t_bins
%size(B_const2)
for i = 1 : num_t_bins,
    B_const(1, i) = B_full(:, i)' * cov_inv * B_full(:, i);
end
%all(all(B_const == B_const2))
%asdnasdks = asndjkasnd

% Precompute values of kernel
h_mat = B_full';
%h_mat = zeros(length(t_sampled), L);
%for i = 1 : length(t_sampled)
%    for n = 1 : L
%        h_mat(i, n) = h(n*T_s - t_sampled(i), args);
%end

setup_obj = struct();
setup_obj.cov_inv = cov_inv;
setup_obj.mu = mu;
setup_obj.t_sampled = t_sampled;
setup_obj.B_full = B_full;
setup_obj.B_const = B_const;
setup_obj.h_mat = h_mat;

end