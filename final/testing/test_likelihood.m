function [t_vec c_vec sig_vec like_vec] = test_likelihood(y, T, args, prev_data)
% Plots best c, sigma, likelihood as a function of spike position t
% Called by run_generic_handle.m

ticID = tic;  % TIMING

K = args.K;
delta_t = args.delta_t;
numiter = args.numiter;
cov_inv = args.setup_obj.cov_inv;
mu = args.setup_obj.mu;
t_sampled = args.setup_obj.t_sampled;
B_full = args.setup_obj.B_full;
B_const = args.setup_obj.B_const;
h = args.h;
delta = args.delta;

% Parameters to be modified
K = 1;
numiter = 3;

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

% Correction for spikes in previous interval
try
    t_prev = prev_data.t - T;
    c_prev = prev_data.c;
    for k = 1 : length(t_prev)
        for n = 1 : L
            y(n) = y(n) - c_prev(k) * h(n*T_s - t_prev(k), args);  % Correction to y for previous spikes
        end
    end
end

% START OF NEW CODE

% t values to plot
t_vec = linspace(0, T, 1000);
c_vec = zeros(size(t_vec));
sig_vec = zeros(size(t_vec));
like_vec = zeros(size(t_vec));

for i = 1 : length(t_vec)

    t = t_vec(i);

    for iter = 1 : numiter

        % Update amplitudes c
        % Maximum likelihood estimate
        % Updates all amplitudes at once

        %dzdc = zeros(K, L);  % dzdc(r, n) = partial derivative dz_n/dc_r  (preallocated)
        for r = 1 : K,
            for n = 1 : L,
                dzdc(r, n) = h(n*T_s - t(r), args);
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
                    d = d + a(r, n) * h(n*T_s - t(k), args);
                end
                A(r, k) = d;
            end
        end

        b = a * y;

        % Apply correction to b arising from non-zero-mean noise in samples
        b = b - a * mu * sigma;

        c = pinv(A) * b;
        
        % Don't allow negative amplitudes
        %c = max(c, 0);

        % Update sigma
        % MLE

        % Compute z
        for n = 1 : L,
            d = 0;
            for k = 1 : K,
                d = d + c(k) * h(n*T_s - t(k), args);
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

    c_vec(i) = c;
    sig_vec(i) = sigma;

    % Compute likelihood

    % Compute z
    for n = 1 : L,
        d = 0;
        for k = 1 : K,
            d = d + c(k) * h(n*T_s - t(k), args);
        end
        z(n) = d;
    end

    % Adjust for sigma
    z = z + sigma / sqrt(2*pi);  % TODO: account for I

    like_vec(i) = exp(-1 / (2 * sigma^2 * (1/2-1/(2*pi)) * delta_t) * (y-z)' * cov_inv * (y-z));

end


return;  % --- BELOW IS OLD CODE (does not run)


for iter = 1:numiter

    %id = tic;

    % Update amplitudes c
    % Maximum likelihood estimate
    % Updates all amplitudes at once

    %dzdc = zeros(K, L);  % dzdc(r, n) = partial derivative dz_n/dc_r  (preallocated)
    for r = 1 : K,
        for n = 1 : L,
            dzdc(r, n) = h(n*T_s - t(r), args);
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
                d = d + a(r, n) * h(n*T_s - t(k), args);
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

    %toc(id)
    %id = tic;

    % Update spike times t_k
    % Finds MLE value by approximating distribution discretely
    % Update times one by one

    for k = 1 : K,  % Estimate t_k

        alpha = y - sigma*mu;
        for kprime = 1 : K,
            if kprime ~= k,
                for n = 1 : L,
                    alpha(n) = alpha(n) - c(kprime) * h(n*T_s - t(kprime), args);
                end
            end
        end

        prob_t = -2*alpha'*cov_inv*B_full + c(k)*B_const;

        % TESTING: impose delta between spikes
        dt_t = t_sampled(2) - t_sampled(1);
        t_samp_temp = t_sampled;
        for kprime = 1 : K
            if kprime ~= k
                rem_ind = find_fast(t_samp_temp, t(kprime) - delta, t(kprime) + delta, dt_t);
                t_samp_temp(rem_ind) = [];
                prob_t(rem_ind) = [];
            end
        end

        [val index] = min(prob_t);
        t(k) = t_samp_temp(index);

    end

    %toc(id)
    %id = tic;

    % Update sigma
    % MLE

    % Compute z
    for n = 1 : L,
        d = 0;
        for k = 1 : K,
            d = d + c(k) * h(n*T_s - t(k), args);
        end
        z(n) = d;
    end

    lambda_0 = (y-z)' * cov_inv * (y-z);
    lambda_1 = 2 * (y-z)' * cov_inv * mu;
    d = (2 * delta_t * (1/2 - 1/(2*pi)) * T);
    alpha_0 = lambda_0 / d;
    alpha_1 = lambda_1 / d;
    sigma = 1/2 * (-alpha_1/L + sqrt(alpha_1^2/L^2 + 8*alpha_0/L));

    %toc(id)

    % TESTING: save intermediate values
    t_inter = [t_inter t];
    c_inter = [c_inter c];
    sig_inter = [sig_inter sigma];

end

[t, c] = sort_t_c(t, c);

elapsed_time = toc(ticID);  % TIMING

% TESTING: plot convergence of params
%{
figure;

subplot(1, 3, 1);
plot(t_inter');
title('t');

subplot(1, 3, 2);
plot(c_inter');
title('c');

subplot(1, 3, 3);
plot(sig_inter);
title('sigma');
%}

end
