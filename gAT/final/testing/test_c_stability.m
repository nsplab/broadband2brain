function [] = test_c_stability(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Tests stability of reconstructed c as sample times change relative to spike

% Params to change
args.K = 1;

max_shift = 0.005;
shift_size = 1;  % In units of dt
shift_vec = 0 : shift_size * get_dt() : max_shift;

t_vec = zeros(size(shift_vec));
c_vec = zeros(size(shift_vec));
sig_vec = zeros(size(shift_vec));

for i = 1 : length(shift_vec)

    % Take samples
    sample_func = eval(['@take_samples_' args.method_name]);
    [y, sample_times] = sample_func(time, dat, t1, T, args);
    %[y] = sample_func(time, dat, t1, T, args);

    % Run reconstruction
    rec_func = eval(['@reconstruct_' args.method_name]);
    [t, c, sigma, elapsed_time] = rec_func(y, T, args, prev_data);

    t_vec(i) = t + shift_vec(i);
    c_vec(i) = c;
    sig_vec(i) = sigma;

    % Shift dat by shift_size indices
    dat = [dat(1 + shift_size : end) zeros(1, shift_size)];

end

figure;
hold on;
plot(shift_vec, t_vec, 'b.-');
plot(shift_vec, c_vec, 'k.-');
plot(shift_vec, sig_vec / 1000, 'm.-');
legend('t', 'c', 'sigma');
xlabel('shift (sec)');
axis tight;
y_lim = ylim;
%ylim([0 y_lim(2)]);
ylim([0 0.04]);
