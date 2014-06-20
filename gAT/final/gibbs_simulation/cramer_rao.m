% Compute Cramer Rao bound on covariance matrix of {c_k, t_k}

% Input:
% n = points at which samples are taken
% t = spike times (true signal)
% c = spike amplitudes (true signal)
% sigmae = std dev of noise
% sigmah = param for gausskernel.m (width of gaussian kernel)

% sampling kernel phi is assumed to be gausskernel.m
% also need derivative of sampling kernel: gaussderiv.m
% see [10] Blu, "Sparse Sampling of Signal Innovations" pg 13

function M = cramer_rao(n, t, c, sigmae, sigmah)

N = length(n);
K = length(t);
R = diag(sigmae^2*ones(1,N), 0);
phi = zeros(N, 2*K);

% left half of phi
for i = 1 : N
    for k = 1 : K
        phi(i, k) = gausskernel(n(i) - t(k), sigmah);
    end
end

% right half of phi
for i = 1 : N
    for k = 1 : K
        phi(i, k + K) = -c(k)*gaussderiv(n(i) - t(k), sigmah);
    end
end

M = (phi'*R^(-1)*phi)^(-1);

end