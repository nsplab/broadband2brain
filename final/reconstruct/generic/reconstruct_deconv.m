function [t, c, sigma, elapsed_time] = reconstruct_deconv(y, T, args, prev_data)
% Fast nonnegative deconvolution using interior point methods
% From Vogelstein paper

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
% sigma - estimate of sigma for gaussian noise
% elapsed_time - how long this function took to run (seconds)

debug = 1;

ticID = tic;

V.dt = T / args.L;
V.est_sig = 1;
V.est_lam = 1;
V.est_gam = 1;
V.est_b = 1;
V.est_a = 1;
V.fast_iter_max = 10;  % iterations for param est
t_y = linspace(V.dt, T, args.L);  % Times at which samples taken
[n_best P_best V C] = fast_oopsi(y, V);

elapsed_time = toc(ticID);

t = zeros(args.K, 1);
c = zeros(args.K, 1);
for i = 1 : args.K
    [val ind] = max(n_best);
    t(i) = t_y(ind);
    c(i) = n_best(ind);
    n_best(ind) = 0;
end

sigma = 0;

if debug
   figure;
   hold on;
   plot(y, 'c');
   stem(10*n_best, 'k');
   plot(C, 'r');
   legend('F', 'n', 'C');
end

end
