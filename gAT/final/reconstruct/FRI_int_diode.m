function [t, c, sigma, elapsed_time] = FRI_int_diode(y, T, K, delta_t, numiter, cov_inv, mu, t_sampled, B_full, B_const)
% Reconstructs spike times and amplitudes using successive integral samples
% Eliminates noise opposite to spikes (assumes spikes negative in x_real)
%   * This can be implemented in analog using a diode
% Takes an approach similar to Gibb's sampling in that it iteratively estimates parameters
% Does not use randomness but instead successively finds MLE for parameters c, t, sigma
% Some large matrices (B_full and B_const) need to be precomputed ahead of time for efficiency
%   * If running in real time make sure they aren't computed every time this is called

% Author: Alex Wein
% Date: Aug, 2011

% Inputs:
% time - time vector for x_real
% x_real - real signal
% L - number of samples
% K - number of spikes
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

ticID = tic;  % TIMING

L = length(y);

% Initial conditions
c = zeros(K, 1);
t = T/2 * ones(K, 1);  %zeros(K, 1);
sigma = 0;

% Preallocated memory
z = zeros(L, 1);
dzdc = zeros(K, L);
alpha = zeros(L, 1);

for iter = 1:numiter

    % Update amplitudes c
    % Maximum likelihood estimate
    % Updates all amplitudes at once

    %dzdc = zeros(K, L);  % dzdc(r, l+1) = partial derivative dz_l/dc_r  (preallocated)
    for r = 1 : K,
        for l = 0 : L-1,
            dzdc(r, l+1) = 1/factorial(l) * (T - t(r))^l;
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
            for l = 0 : L-1,
                d = d + a(r, l+1)/factorial(l) * (T-t(k))^l;
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
    %c

    % Update spike times t_k
    % Finds MLE value by approximating distribution discretely
    % Update times one by one

    for k = 1 : K,  % Estimate t_k

        alpha = y - sigma*mu;
        for kprime = 1 : K,
            if kprime ~= k,
                for l = 0 : L-1,
                    alpha(l+1) = alpha(l+1) - 1/factorial(l) * c(kprime) * (T - t(kprime))^l;
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
    %t

    % Update sigma
    % MLE

    % Compute z
    for l = 0 : L-1,
        d = 0;
        for k = 1 : K,
            d = d + c(k) * (T - t(k))^l;
        end
        z(l+1) = d / factorial(l);
    end

    lambda_0 = (y-z)' * cov_inv * (y-z);
    lambda_1 = 2 * (y-z)' * cov_inv * mu;
    d = (2 * delta_t * (1/2 - 1/(2*pi)) * T);
    alpha_0 = lambda_0 / d;
    alpha_1 = lambda_1 / d;
    sigma = 1/2 * (-alpha_1/L + sqrt(alpha_1^2/L^2 + 8*alpha_0/L));
    %sigma

end

[t, c] = sort_t_c(t, c);

elapsed_time = toc(ticID);  % TIMING

end
