function [t, c, sigma, elapsed_time] = reconstruct_int(y, T, args, prev_data)
% Annihilating filter
% L = 2*K+1

K = args.K;

ticID = tic;

% Form matrix S
S = zeros(K+1, K+1);
for i = 1 : K+1,
    for j = 1 : K+1,
        S(i, j) = y(K + i - j + 1);
    end
end

% Singular value decomposition
[U, Sigma, V] = svd(S);

% Extract annihilating filter coefficients
a = V(:,K+1);

%a

% Find roots of A(z), the z-transform of the annihilating filter
% Note that multiplying by the correct number of factors of z turns A(z)
% into a polynomial with coefficients a_0, a_1, ... , a_K. Therefore we can
% use matlab's polynomial root function.
u = roots(a);
j = 1; % Number of (reasonable) spikes recovered is j-1
t = []; % Reconstructed times
u_ind = []; % u_ind[i] = index in u associated with t(i)
for i = 1 : length(u) %K
    % Only keep roots that are reasonable (i.e. real and between 0 and T)
    if imag(u(i)) == 0 && u(i) > 0 && u(i) < T,
        t(j) = T - u(i);
        u_ind(j) = i;
        j = j + 1;
    end
end

%u

% Reconstruct spike amplitudes
% Solve U*c = s (step 3 of the algorithm in the Srinivasan paper)
% TODO: what if not all roots found? (then u still has length K but some values are bad)
% Construct U from u
U = zeros(K, K);
U(1, 1:K) = ones(1, 1:K);
for i = 2 : K,
    U(i, 1:K) = U(i-1, 1:K) .* u';
end
c = real(pinv(U) * y(1:K, 1)); % Amplitudes

elapsed_time = toc(ticID);  % TIMING

sigma = 0;
t = sort(t)';

end
