function [res next_data] = run_sinc_diode_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Applies diode to data
% Takes sinc samples
% Runs reconstruction method and returns spike times/amplitudes
% res(1, :) = times
% res(2, :) = amplitudes
% res(3, 1) = sigma
% Should be run with small window around data (see iterate_itervals.m) and with data that is already highpassed
% args is a struct with fields L, K, B, delta_t, numiter, cov_inv, mu, t_sampled, B_full, B_const

%ticID_sinc = tic;

% Apply ideal diode
dat = ideal_diode(dat);

% TESTING: fake data
%dat = dat * 0;
%dat(round(length(dat)/2)) = 500;

% Take sinc samples
%[y, sample_times] = sample_sinc(time, dat, t1, T, args.B, args.L);
[y, sample_times] = convolve_fast(dat, t1);  % Optimized

%sinc_time_mid = toc(ticID_sinc)

% Run reconstruction
[t, c, sigma, elapsed_time] = FRI_sinc_diode(y, T, args.K, args.B, args.delta_t, args.numiter, args.cov_inv, args.mu, args.t_sampled, args.B_full, args.B_const);
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

% TESTING: plot
%{
sinc_t = time(1) : 0.0001 : time(end);
sinc_x = zeros(size(sinc_t));
for i = 1 : length(sinc_x)
    sinc_x(i) = normalized_sinc(sinc_t(i) - mean(sinc_t), [args.B]);
end
figure;
hold on;
plot(time, dat, 'c');
plot(sample_times, y*5, 'm.-');
stem(t, c*10000, 'k');
plot([t1, t1+T], [-5, -5], 'g.');
plot(sinc_t, sinc_x, 'r');
legend('data', 'samples', 'reconstructed', 'interval', 'sinc');
xlabel('time (sec)');
title(['sigma = ' num2str(sigma)]);
%}

%sinc_time = toc(ticID_sinc)

next_data = 0;

end
