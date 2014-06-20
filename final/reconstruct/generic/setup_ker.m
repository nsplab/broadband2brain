function [setup_obj] = setup_ker(T, window, args)
% Precomputes values for generic kernel
% args.h is function handler for kernel
% args.h_args is args for kernel

L = args.L;
h = args.h;
t_resolution = args.t_resolution;

% Set up for convolution
convolve_setup(T, window, L, h, args);

T_s = T / L;  % Spacing of samples

% Covariance matrix for samples y
% Samples assumed to follow a multivariate gaussian distribution
cov = zeros(L, L);
for i = 1 : L,
    for j = 1 : L,
        cov(i, j) = get_cov(abs(i-j), T_s, h, args);  % Was wrong, T should be T_s
    end
end

cov_inv = cov^(-1);

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

function cov = get_cov(delta_n, T_s, h, args)
% NOTE: this needs to be updated for each new kernel
% Gives value of integral for covariance term

if strcmp(func2str(h), 'normalized_sinc')
    cov = h(delta_n * T_s, args);
elseif strcmp(func2str(h), 'RC_lowpass')
    cov = h(delta_n * T_s, args) / 2;
elseif strcmp(func2str(h), 'RLC')
    s1 = args.s1;
    s2 = args.s2;
    cov = (s1*s2/(s2-s1))^2 * (exp(s1*delta_n*T_s) * (1/(s1+s2)-1/(2*s1)) + exp(s2*delta_n*T_s) * (1/(s1+s2)-1/(2*s2)));
end

end
