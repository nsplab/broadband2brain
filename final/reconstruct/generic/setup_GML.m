function [setup_obj] = setup_GML(T, window, args)
% Precomputes values for generic kernel
% args.h is function handler for kernel
% args.h_args is args for kernel

L = args.L;
h = args.h;
J = args.J;

ns = linspace(0, T, L);  % points at which samples are taken

% Covariance matrix for samples y
% Samples assumed to follow a multivariate gaussian distribution
cov_inv = eye(L);

% Compute mu, used to correct for non-zero-mean noise in samples
mu = zeros(L, 1);

% Discrete values of t
t_sampled = linspace(0, T, J);

% Precompute values for grid search to estimate t_k
B_full = zeros(L, J);
for n = 1 : L,
    for i = 1 : J,
        B_full(n, i) = h(ns(n) - t_sampled(i), args);
    end
end
B_const = zeros(1, J);
for i = 1 : J,
    B_const(1, i) = B_full(:, i)' * cov_inv * B_full(:, i);
end

% Precompute values of kernel
% THIS IS REDUNDANT!!! Same as B_full
%{
h_mat = zeros(length(t_sampled), L);
for i = 1 : length(t_sampled)
    for n = 1 : L
        h_mat(i, n) = h(ns(n) - t_sampled(i), args);
    end
end
%}

setup_obj = struct();
setup_obj.cov_inv = cov_inv;
setup_obj.mu = mu;
setup_obj.t_sampled = t_sampled;
setup_obj.B_full = B_full;
setup_obj.B_const = B_const;
%setup_obj.h_mat = h_mat;

end
