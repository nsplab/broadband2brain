% Compute Cramer Rao bound on covariance matrix of {c_k, t_k}

% Input:
% n = points at which samples are taken
% t = spike times (true signal)
% c = spike amplitudes (true signal)
% sigmae = std dev of noise
% sigmah = param for gausskernel.m (width of gaussian kernel)

% sampling kernel phi is assumed to be sinc_kernel.m (periodic sinc)
% also need derivative of sampling kernel: sinc_deriv.m
% see [10] Blu, "Sparse Sampling of Signal Innovations" pg 13

function M = cramer_rao_psinc(n, t, c, sigmae, B, tau)

N = length(n);
K = length(t);
R = diag(sigmae^2*ones(1,N), 0);
phi = zeros(N, 2*K);

% left half of phi
for i = 1 : N
    for k = 1 : K
        phi(i, k) = sinc_kernel(n(i) - t(k), B, tau);
    end
end

% right half of phi
for i = 1 : N
    for k = 1 : K
        phi(i, k + K) = -c(k)*sinc_deriv(n(i) - t(k), B, tau);
    end
end

M = (phi'*R^(-1)*phi)^(-1);

end