function [res next_data] = run_int_diode_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Applies diode to data
% Takes successive integral samples
% Runs reconstruction method and returns spike times/amplitudes
% res(1, :) = times
% res(2, :) = amplitudes
% res(3, 1) = sigma
% Should be run with no window around data (see iterate_itervals.m) and with data that is already highpassed
% args is a struct with fields L, K, delta_t, numiter, cov_inv, mu, t_sampled, B_full, B_const

% Apply ideal diode
dat = ideal_diode(dat);

% Take successive integral samples
y = succ_int(dat, get_dt(), args.L);

% Run reconstruction
[t, c, sigma, elapsed_time] = FRI_int_diode(y, T, args.K, args.delta_t, args.numiter, args.cov_inv, args.mu, args.t_sampled, args.B_full, args.B_const);
t = t + t1;

% Limit spike times to [seg_start, seg_end]
ind = find(t >= seg_start & t < seg_end);
t = t(ind, 1);
c = c(ind, 1);

% Output
res = zeros(3, length(t));

res(1, :) = t';
res(2, :) = c';
res(3, :) = sigma * ones(size(res(3, :)));

next_data = 0;

end
