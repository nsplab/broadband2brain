function [t, c, m_w, elapsed_time] = reconstruct_ML(y, T, args, prev_data)
% Reconstructs spike times and amplitudes using samples taken by kernel
% Iteratively estimates parameters [c, t, sigma] using MLE
% Spikes must have positive amplitude

% Author: Alex Wein
% Date: Aug, 2011

% Inputs:
% y - samples taken by convolution with sinc function (see sinc_sample.m)
% T - length of time interval
% args
% prev_data - data from previous interval, contains previous spike times and amplitudes (prev_data.t, prev_data.c)
% --- precomputed values ---
% cov_inv
% mu
% t_sampled
% B_full
% B_const
% h_mat

% Outputs:
% t - spike times
% c - spike amplitudes
% m_w - estimate for mean of noise (after rectification)
% elapsed_time - how long this function took to run (seconds)

ticID = tic;  % TIMING

debug = 0;  % Makes plots if set to 1

K = args.K;
delta_t = args.delta_t;
numiter = args.numiter;
cov_inv = args.setup_obj.cov_inv;
mu = args.setup_obj.mu;
t_sampled = args.setup_obj.t_sampled;
B_full = args.setup_obj.B_full;
B_const = args.setup_obj.B_const;
h_mat = args.setup_obj.h_mat;
h = args.h;
delta = args.delta;

L = length(y);
T_s = T / L;  % Spacing of samples

% Initial conditions
c = zeros(K, 1);
t_ind = round(length(t_sampled)/2) * ones(K, 1);
t = t_sampled(t_ind)';
m_w = 0;

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

% TESTING: save intermediate values
%t_inter = [];
%c_inter = [];
%sig_inter = [];

for iter = 1:numiter

    %disp('start iter');
    %id = tic;

    % Update amplitudes c
    % Maximum likelihood estimate
    % Updates all amplitudes at once

    %dzdc = zeros(K, L);  % dzdc(r, n) = partial derivative dz_n/dc_r  (preallocated)
    %for r = 1 : K,
    %    for n = 1 : L,
    %        dzdc(r, n) = h(n*T_s - t(r), args);  % TODO: optimize, precompute
    %    end
    %end
    dzdc = h_mat(t_ind, :);

    %a = zeros(K, L);  % a_rl
    %for r = 1 : K,
    %    a(r, :) = dzdc(r, :) * cov_inv;
    %end
    a = dzdc * cov_inv;

    A = zeros(K, K);  % A_rk
    for r = 1 : K,
        for k = 1 : K,
            d = 0;
            for n = 1 : L,
                d = d + a(r, n) * h_mat(t_ind(k), n);
            end
            A(r, k) = d;
        end
    end

    b = a * y;

    % Apply correction to b arising from non-zero-mean noise in samples
    b = b - a * mu * m_w;

    c = pinv(A) * b;
    
    % Don't allow negative amplitudes
    c = max(c, 0);

    % TESTING: make plots to analyze this step
    %{
    if debug == 1
        if K == 2
            this_k = 2;
            sig_step = 0.001;
            num_plot = 20;
            sig_vec = c(this_k) - sig_step * num_plot : sig_step : c(this_k) + sig_step * num_plot;
            L_vec = zeros(size(sig_vec));
            for i = 1 : length(sig_vec)
                c_val = sig_vec(i);
                mu_vec = zeros(L, 1);
                for n = 1 : L,
                    d = 0;
                    for k = 1 : K,
                        if k == this_k
                            d = d + c_val * h(n*T_s - t(k), args);
                        else
                            d = d + c(k) * h(n*T_s - t(k), args);
                        end
                    end
                    mu_vec(n) = d;
                end
                mu_vec = mu_vec + m_w;
                tau_stuff(i) = (y-mu_vec)' * cov_inv * (y-mu_vec);
                %exp_arg(i) = -1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec));
                %L_vec(i) = 1/sigma^L * exp(-1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec)));
            end
            figure;
            hold on;
            %plot(sig_vec, L_vec, 'b');
            plot(sig_vec, tau_stuff, 'g');
            %plot(sig_vec, exp_arg, 'y');
            plot([c(this_k) c(this_k)], ylim, 'r');
            xlabel('c_k');
            title('c');
            legend('L (likelihood)', 'c (chosen value)');
        end
    end
    %}
    % END TESTING

    %toc(id)
    %id = tic;

    % Update spike times t_k
    % Finds MLE value by approximating distribution discretely
    % Update times one by one

    for k = 1 : K,  % Estimate t_k

        alpha = y - m_w;
        for kprime = 1 : K,
            if kprime ~= k,
                for n = 1 : L,
                    alpha(n) = alpha(n) - c(kprime) * h_mat(t_ind(kprime), n);
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

        % Impose delta between spikes
        dt_t = t_sampled(2) - t_sampled(1);
        t_samp_temp_ind = 1 : length(t_sampled);
        t_samp_temp = t_sampled;
        for kprime = 1 : K
            if kprime ~= k
                % TESTING *** remove impose delta: rem_ind = [];
                rem_ind = find_fast(t_samp_temp, t(kprime) - delta, t(kprime) + delta, dt_t);
                t_samp_temp_ind(rem_ind) = [];
                t_samp_temp(rem_ind) = [];
                prob_t(rem_ind) = [];
            end
        end

        [val index] = min(prob_t);
        t_ind(k) = t_samp_temp_ind(index);
        t(k) = t_samp_temp(index);
        
        % TESTING: make plots to analyze this step
        %{
        if debug == 1
            sig_step = 0.00001;
            num_plot = 200;
            sig_vec = t(k) - sig_step * num_plot : sig_step : t(k) + sig_step * num_plot;
            L_vec = zeros(size(sig_vec));
            for i = 1 : length(sig_vec)
                t_val = sig_vec(i);
                mu_vec = zeros(L, 1);
                for n = 1 : L,
                    d = 0;
                    for k_prime = 1 : K,
                        if k == k_prime
                            d = d + c(k) * h(n*T_s - t_val, args);
                        else
                            d = d + c(k) * h(n*T_s - t(k_prime), args);
                        end
                    end
                    mu_vec(n) = d;
                end
                mu_vec = mu_vec + m_w;
                tau_stuff(i) = (y-mu_vec)' * cov_inv * (y-mu_vec);
                %exp_arg(i) = -1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec));
                %L_vec(i) = 1/sigma^L * exp(-1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec)));
            end
            figure;
            hold on;
            %plot(sig_vec, L_vec, 'b');
            plot(sig_vec, tau_stuff, 'g');
            %plot(sig_vec, exp_arg, 'y');
            plot([t(k), t(k)], ylim, 'r');
            xlabel('t_k');
        end
        %}
        % END TESTING
    
        % C ITER
        %{
        dzdc = h_mat(t_ind, :);

        %a = zeros(K, L);  % a_rl
        %for r = 1 : K,
        %    a(r, :) = dzdc(r, :) * cov_inv;
        %end
        a = dzdc * cov_inv;

        A = zeros(K, K);  % A_rk
        for r = 1 : K,
            for k = 1 : K,
                d = 0;
                for n = 1 : L,
                    d = d + a(r, n) * h_mat(t_ind(k), n);
                end
                A(r, k) = d;
            end
        end

        b = a * y;

        % Apply correction to b arising from non-zero-mean noise in samples
        b = b - a * mu * m_w;

        c = pinv(A) * b;

        % Don't allow negative amplitudes
        c = max(c, 0);
        %}
        % END C ITER

    end

    %toc(id)
    %id = tic;

    % Update m_w
    % MLE

    d = 0;
    for n = 1 : L,
        for k = 1 : K,
            d = d + c(k) * h(n*T_s - t(k), args);
        end
    end
    
    m_w = (mu' * cov_inv * y - d) / L;

    %lambda_0 = (y-z)' * cov_inv * (y-z);
    %lambda_1 = 2 * (y-z)' * cov_inv * mu;
    %d = (2 * delta_t * (1/2 - 1/(2*pi)));  % removed factor of T (typo)
    %alpha_0 = lambda_0 / d;
    %alpha_1 = lambda_1 / d;
    %sigma = 1/2 * (-alpha_1/L + sqrt(alpha_1^2/L^2 + 8*alpha_0/L));
    
    % TESTING: make plots to analyze this step
    %{
    if debug == 1
        if K == 2
            sig_step = 0.5;
            num_plot = 20;
            sig_vec = m_w - sig_step * num_plot : sig_step : m_w + sig_step * num_plot;
            L_vec = zeros(size(sig_vec));
            mu_vec = zeros(L, 1);    
            tau_stuff = zeros(size(sig_vec));
            for i = 1 : length(sig_vec)
                for n = 1 : L,
                    d = 0;
                    for k = 1 : K,
                        d = d + c(k) * h(n*T_s - t(k), args);
                    end
                    mu_vec(n) = d;
                end
                s = sig_vec(i);
                mu_vec = mu_vec + s;
                tau_stuff(i) = (y-mu_vec)' * cov_inv * (y-mu_vec);
                %exp_arg(i) = -1/(2*s^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec));
                %L_vec(i) = 1/s^L * exp(-1/(2*s^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec)));
            end
            figure;
            hold on;
            %plot(sig_vec, L_vec, 'b');
            plot(sig_vec, tau_stuff, 'g');
            %plot(sig_vec, exp_arg / 100, 'y');
            plot([m_w m_w], ylim, 'r');
            title('m_w');
            legend('L (likelihood)', 'm_w (chosen value)');
        end
    end
    %}
    % END TESTING
    
    %toc(id)

    % TESTING: save intermediate values
    %t_inter = [t_inter t];
    %c_inter = [c_inter c];
    %sig_inter = [sig_inter sigma];
    
    %%% ***** Extra c iteration ******
    
    %{
    
    % Update amplitudes c
    % Maximum likelihood estimate
    % Updates all amplitudes at once

    %dzdc = zeros(K, L);  % dzdc(r, n) = partial derivative dz_n/dc_r  (preallocated)
    %for r = 1 : K,
    %    for n = 1 : L,
    %        dzdc(r, n) = h(n*T_s - t(r), args);  % TODO: optimize, precompute
    %    end
    %end
    dzdc = h_mat(t_ind, :);

    %a = zeros(K, L);  % a_rl
    %for r = 1 : K,
    %    a(r, :) = dzdc(r, :) * cov_inv;
    %end
    a = dzdc * cov_inv;

    A = zeros(K, K);  % A_rk
    for r = 1 : K,
        for k = 1 : K,
            d = 0;
            for n = 1 : L,
                d = d + a(r, n) * h_mat(t_ind(k), n);
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
    
    %}
    % End extra c iter
    
    % TESTING: make plots to analyze this step
    %{
    if K == 2 && iter == numiter
        this_k = 2;
        sig_step = 0.001;
        num_plot = 20;
        sig_vec = c(this_k) - sig_step * num_plot : sig_step : c(this_k) + sig_step * num_plot;
        L_vec = zeros(size(sig_vec));
        for i = 1 : length(sig_vec)
            c_val = sig_vec(i);
            mu_vec = zeros(L, 1);
            for n = 1 : L,
                d = 0;
                for k = 1 : K,
                    if k == this_k
                        d = d + c_val * h(n*T_s - t(k), args);
                    else
                        d = d + c(k) * h(n*T_s - t(k), args);
                    end
                end
                mu_vec(n) = d;
            end
            mu_vec = mu_vec + sigma / sqrt(2 * pi);
            tau_stuff(i) = (y-mu_vec)' * cov_inv * (y-mu_vec);
            norm(i) = (y-mu_vec)' * (y-mu_vec);
            exp_arg(i) = -1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec));
            L_vec(i) = 1/sigma^L * exp(-1/(2*sigma^2 * (1/2-1/(2*pi)) * delta_t) * ((y-mu_vec)' * cov_inv * (y-mu_vec)));
        end
        figure;
        hold on;
        %plot(sig_vec, L_vec, 'b');
        plot(sig_vec, tau_stuff, 'g');
        plot(sig_vec, norm, 'm');
        %plot(sig_vec, exp_arg, 'y');
        plot([c(this_k) c(this_k)], ylim, 'r');
        xlabel('c_k');
    end
    %}
    % END TESTING

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
