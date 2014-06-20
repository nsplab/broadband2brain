function [t, c, sigma, elapsed_time] = FRI_sinc_diode(y, T, K, B, delta_t, numiter, cov_inv, mu, t_sampled, B_full, B_const)
% Reconstructs spike times and amplitudes using samples taken by sinc kernel with diode
% Iteratively estimates parameters [c, t, sigma] using MLE
% Samples should be taken with sinc_sample.m after eliminating negative noise with diode
% Spikes must have positive amplitude

% Author: Alex Wein
% Date: Aug, 2011

% Inputs:
% y - samples taken by convolution with sinc function (see sinc_sample.m)
% T - length of time interval
% K - number of spikes
% B - bandwidth of sinc
% delta_t - noise modeled as i.i.d. gaussian with this bin width, should be set to about dt
% numiter - number of iterations
% --- precomputed values ---
% cov_inv
% mu
% t_sampled
% B_full
% B_const

% Outputs:
% t - spike times
% c - spike amplitudes
% sigma - estimate of sigma for gaussian noise
% elapsed_time - how long this function took to run (seconds)

ticID = tic;  % TIMING

L = length(y);
T_s = T / L;  % Spacing of samples

% Initial conditions
c = zeros(K, 1);
t = T/2 * ones(K, 1);
sigma = 0;

% Preallocated memory
z = zeros(L, 1);
dzdc = zeros(K, L);
alpha = zeros(L, 1);

for iter = 1:numiter

    % Update amplitudes c
    % Maximum likelihood estimate
    % Updates all amplitudes at once

    %dzdc = zeros(K, L);  % dzdc(r, n) = partial derivative dz_n/dc_r  (preallocated)
    for r = 1 : K,
        for n = 1 : L,
            dzdc(r, n) = normalized_sinc((n-1/2)*T_s - t(r), B);
        end
    end

    a = zeros(K, L);  % a_rl
    for r = 1 : K,
        a(r, :) = dzdc(r, :) * cov_inv;
    end

    A = zeros(K, K);  % A_rk
    for r = 1 : K,
        for k = 1 : K,
            d = 0;
            for n = 1 : L,
                d = d + a(r, n) * normalized_sinc((n-1/2)*T_s - t(k), B);
            end
            A(r, k) = d;
        end
    end

    b = a * y;

    % Apply correction to b arising from non-zero-mean noise in samples
    b = b - a * mu * sigma;

    c = pinv(A) * b;
    
    % Don't allow negative amplitudes
    c = max(c, 0);

    % Update spike times t_k
    % Finds MLE value by approximating distribution discretely
    % Update times one by one

    for k = 1 : K,  % Estimate t_k

        alpha = y - sigma*mu;
        for kprime = 1 : K,
            if kprime ~= k,
                for n = 1 : L,
                    alpha(n) = alpha(n) - c(kprime) * normalized_sinc((n-1/2)*T_s - t(kprime), B);
                end
            end
        end

        prob_t = -2*alpha'*cov_inv*B_full + c(k)*B_const;

        % TESTING: plot distribution
        %{
        figure;
        plot(t_sampled, prob_t, 'b.-')
        title(['iter ' int2str(iter) '  k = ' int2str(k)]);
        %}

        [val index] = min(prob_t);
        t(k) = t_sampled(index);

    end

    % Update sigma
    % MLE

    % Compute z
    for n = 1 : L,
        d = 0;
        for k = 1 : K,
            d = d + c(k) * normalized_sinc((n-1/2)*T_s - t(k), B);
        end
        z(n) = d;
    end

    lambda_0 = (y-z)' * cov_inv * (y-z);
    lambda_1 = 2 * (y-z)' * cov_inv * mu;
    d = (2 * delta_t * (1/2 - 1/(2*pi)) * T);
    alpha_0 = lambda_0 / d;
    alpha_1 = lambda_1 / d;
    sigma = 1/2 * (-alpha_1/L + sqrt(alpha_1^2/L^2 + 8*alpha_0/L));

end

[t, c] = sort_t_c(t, c);

elapsed_time = toc(ticID);  % TIMING

end
