% Runs run_FRI_int_diode.m
% Plots spike train

close all;
clear all;

ticID = tic;

% Parameters
data_set = 1;
electrodes = [9]; %channels_to_use(data_set);

[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = data_start;
end_time = start_time + 10;
cond = 'ih';
T = 0.02;
L = 3;
K = 1;
delta_t = get_dt();
numiter = 3;
t_resolution = T / 100;
tolerance = 0.01;

% Estimate optimal thresholds   TODO: only necessary if parameters have changed, can be commented otherwise
%est_rej_thres(data_set, electrodes, cond, cond, T, L, K, delta_t, numiter, t_resolution, tolerance);

for electrode = electrodes

    % FRI_int_diode spikes
    [t, c, t1, sigma] = run_FRI_int_diode(data_set, [electrode], start_time, end_time, cond, start_time, T, L, K, delta_t, numiter, t_resolution);
    t = t{1};
    c = c{1};
    t1 = [t1 t1(end)+T];
    sigma = sigma{1};

    % 2 methods for getting real spikes: 'real' and threshold

    % Real spikes
    %spikes = real_spikes(data_set, electrode, start_time, end_time, 1);

    % Analog threshold spikes
    spikes = run_threshold(data_set, [electrode], start_time, end_time, cond);
    spikes = spikes{1};

    % Get optimal threshold
    thres = choose_threshold(data_set, electrode, cond, 1);  % Must load, compute using est_rej_thres.m

    % Apply threshold to reject reconstructed spikes
    keep_ind = find(c >= thres);
    keep_t = t(keep_ind);
    keep_c = c(keep_ind);
    rej_ind = find(c < thres);
    rej_t = t(rej_ind);
    rej_c = c(rej_ind);

    % Compare spikes
    [tp, fp, fn, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(keep_t, spikes, tolerance);

    scale = 1/(get_dt()*10);

    % Plot results
    fig = figure;
    hold on;
    [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
    plot(time, dat, 'c');  % Real data
    stem(keep_t, keep_c*scale, 'k');  % Reconstructed spikes
    stem(rej_t, rej_c*scale, 'Color', [0.5, 0.5, 0.5]);  % Rejected spikes
    plot(t1, ones(size(t1))*-25, 'b.');  % Sampling times
    plot(spikes, ones(size(spikes))*-15, 'k.');  % Real spikes
    plot([time(1), time(end)], scale*[thres, thres], 'g--');  % Threshold
    stem(t, sigma, 'm');  % Plot sigma

    % Plot false positives, etc.
    plot(keep_t(tp), ones(size(tp))*-55, 'g.');
    plot(keep_t(fp), ones(size(fp))*-45, 'm.');
    plot(spikes(fn), ones(size(fn))*-35, 'r.');

    leg = legend('data', 'reconstructed spikes', 'rejected spikes', 'sampling times', 'real spikes', 'threshold', 'sigma', 'true positives', 'false positives', 'false negatives');
    set(leg, 'Location', 'SouthOutside');
    xlabel('time (sec)');

    % Output error rates
    title({['data set ' int2str(data_set) ', electrode ' int2str(electrode)], ['true positive rate: ' num2str(tp_rate)], ['false positive rate: ' num2str(fp_rate)], ['miss rate: ' num2str(fn_rate)], ['error rate: ' num2str(error_rate)]});

    saveas(fig, ['~/Desktop/spk_cmp_' int2str(data_set) '-' int2str(electrode) '.pdf']);

end

toc(ticID)
