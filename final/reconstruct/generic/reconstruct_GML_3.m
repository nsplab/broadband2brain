function [t, c, sigmae, elapsed_time, t_inter, c_inter, sig_inter] = reconstruct_GML_2(y, T, args, prev_data)

% GML 2
% - update t_k, c
% - c > 0
% - delta

% GML: ML algorithm for comparison with Gibbs sampling on simulated data
% Reconstructs spike times and amplitudes using samples taken by kernel
% Iteratively estimates parameters [c, t, sigmae] using MLE

% Author: Alex Wein
% Date: July, 2012

% Inputs:
% y - samples (gaussian kernel)
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
numiter = args.numiter;
cov_inv = args.setup_obj.cov_inv;
t_sampled = args.setup_obj.t_sampled;
B_full = args.setup_obj.B_full;
B_const = args.setup_obj.B_const;
h_mat = args.setup_obj.h_mat;
h = args.h;
delta = args.delta;

L = length(y);
ns = linspace(0, T, L);  % points at which samples are taken

% Initial conditions
c = zeros(K, 1);
t_ind = round(length(t_sampled)/2) * ones(K, 1);
t = t_sampled(t_ind)';

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

% Save intermediate values
t_inter = [];
c_inter = [];
sig_inter = [];

for iter = 1 : numiter

    % Update spike times t_k
    % Finds MLE value by approximating distribution discretely
    % Update times one by one

    for k = 1 : K,  % Estimate t_k

        alpha = y;% - m_w;
        for kprime = 1 : K,
            if kprime ~= k,
                for n = 1 : L,
                    alpha(n) = alpha(n) - c(kprime) * h_mat(t_ind(kprime), n);
                end
            end
        end

        prob_t = -2*alpha'*cov_inv*B_full + c(k)*B_const;

        % Impose delta between spikes
        dt_t = t_sampled(2) - t_sampled(1);
        t_samp_temp_ind = 1 : length(t_sampled);
        t_samp_temp = t_sampled;
        for kprime = 1 : K
            if kprime ~= k
                % TESTING *** to remove impose delta: rem_ind = [];
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
        if debug == 1
            t_vec = t(k)-1 : 0.01 : t(k)+1;
            L_vec = zeros(size(t_vec));
            for i = 1 : length(L_vec)
                tprime = t;
                tprime(k) = t_vec(i);
                L_vec(i) = lh(tprime, c, 100, y, ns, K, L, args.sigmah);
            end
            figure;
            hold on;
            plot(t_vec, L_vec, 'g');
            plot([t(k) t(k)], ylim, 'r');
            axis tight;
            plot(t_sampled, min(L_vec)*ones(size(t_sampled)), 'b.');
            title('t');
        end
        % END TESTING
        
        % Update amplitudes c
        % Maximum likelihood estimate
        % Updates all amplitudes at once

        dzdc = h_mat(t_ind, :);
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
        c = pinv(A) * b;

        % Don't allow negative amplitudes
        c = max(c, 0);

        % TESTING: make plots to analyze this step
        if debug == 1
            kprime = 2;
            c_vec = c(kprime)-1 : 0.01 : c(kprime)+1;
            L_vec = zeros(size(c_vec));
            for i = 1 : length(L_vec)
                cprime = c;
                cprime(kprime) = c_vec(i);
                L_vec(i) = lh(t, cprime, 100, y, ns, K, L, args.sigmah);
            end
            figure;
            hold on;
            plot(c_vec, L_vec, 'g');
            plot([c(kprime) c(kprime)], ylim, 'r');
            title('c');
        end
        % END TESTING

    end
    
    c_inter = [c_inter c];  % save intermediate results
    t_inter = [t_inter t];  % save intermediate results

end

[t, c] = sort_t_c(t, c);

% Estimate sigmae
z = zeros(L, 1);
for k = 1 : K
    for i = 1 : L
        z(i) = z(i) + c(k) * exp(-(ns(i)-t(k))^2/(2*args.sigmah^2));
    end
end
M = 0;
for i = 1 : L
    M = M + (y(i)-z(i))^2;
end
sigmae = sqrt(M / (L+1));

sig_inter = [sig_inter sigmae];  % save intermediate results

% TESTING: make plots to analyze this step
if debug == 1
    s_vec = sigmae-1 : 0.01 : sigmae+1;
    L_vec = zeros(size(s_vec));
    for i = 1 : length(L_vec)
        L_vec(i) = lh(t, c, s_vec(i), y, ns, K, L, args.sigmah);
    end
    figure;
    hold on;
    plot(s_vec, L_vec, 'g');
    plot([sigmae sigmae], ylim, 'r');
    title('sigma_e');
end
% END TESTING

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

% Calculates data likelihood
% Only used for testing purposes
function L = lh(t, c, sigmae, y, ns, K, L, sigmah)

z = zeros(L, 1);
for k = 1 : K
    for i = 1 : L
        z(i) = z(i) + c(k) * exp(-(ns(i)-t(k))^2/(2*sigmah^2));
    end
end
M = 0;
for i = 1 : L
    M = M + (y(i)-z(i))^2;
end

%L = M;
L = 1/sigmae^(L+1) * exp(-1/(2*sigmae^2)*M);

end
