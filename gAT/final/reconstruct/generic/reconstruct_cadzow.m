% Implementation of [10] Blu "Sparse Sampling of Signal Innovations"
% Annihilating filter method with cadzow denoising
% But adapted for Gaussian kernel using Vetterli "Sampling Signals with
% Finite Rate of Innovation"
% Gaussian kernel: exp(-t^2/(2*sigmah^2))

function [t, c, sigmae, elapsed_time] = reconstruct_cadzow(y, T, args)

N = length(y);
n = linspace(0, T, N);  % sampling points
T_s = n(2)-n(1);  % called T in Vetterli paper
K = args.K;
sigmah = args.sigmah;

ticID = tic;

%y

S_vec = y .* exp((0:N-1)'.^2*T_s^2/(2*sigmah^2));

%S_vec

S_vec = cadzow_denoise(S_vec, K);

%S_vec_noiseles = S_vec

% Form matrix S
S = zeros(N-K, K+1);
for i = 1 : N-K,
    for j = 1 : K+1,
        S(i, j) = S_vec(K + i - j + 1);
    end
end

% Singular value decomposition
[U, Sigma, V] = svd(S);

% Extract annihilating filter coefficients
a = V(:,K+1);

% Find roots of A(z), the z-transform of the annihilating filter
% Note that multiplying by the correct number of factors of z turns A(z)
% into a polynomial with coefficients a_0, a_1, ... , a_K. Therefore we can
% use matlab's polynomial root function.
u = roots(a);
t = log(u)*sigmah^2/T_s;
%{
j = 1; % Number of (reasonable) spikes recovered is j-1
t = []; % Reconstructed times
u_ind = []; % u_ind[i] = index in u associated with t(i)
for i = 1 : length(u) %K
    % Only keep roots that are reasonable (i.e. real and between 0 and T)
    if imag(u(i)) == 0 && u(i) > 0 && u(i) < T,
        t(j) = u(i);
        u_ind(j) = i;
        j = j + 1;
    end
end
%}

% Reconstruct spike amplitudes
% Solve U*c = s (step 3 of the algorithm in the Srinivasan paper)
% TODO: what if not all roots found? (then u still has length K but some values are bad)
% Construct U from u
U = zeros(N, K);
U(1, 1:K) = ones(1, 1:K);
for i = 2 : N,
    U(i, 1:K) = U(i-1, 1:K) .* u';
end
a = real(pinv(U) * S_vec); % Amplitudes
c = a .* exp(t.^2/(2*sigmah^2));

elapsed_time = toc(ticID);

sigmae = 0;
[t, c] = sort_t_c(t', c');

end
