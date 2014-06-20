function [] = test_analog()
% Tests analog FRI

close all;

% Parameters
data_set = 1;
electrodes = channels_to_use(data_set);
cond = 'raw';
duration = 2;

method_name = 'analog'; %'gibbs' for conference
paramID = 1;


[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = data_start;
end_time = start_time + duration;

[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);

tolerance = 0.0025;  % for false-positives etc

for electrode = electrodes

    % Reconstruct spikes
    [t1, result] = run_generic(data_set, [electrode], [start_time], [end_time], hp_handle, diode_handle, method_name, args, start_time, T, window);
    [t, c, t1, sigma] = extract_result_single(t1, result);
    t = t{1};
    c = c{1};
    t1 = [t1 t1(end)+T];
    sigma = sigma{1};

    % Real spikes
    spikes = real_spikes(data_set, electrode, start_time, end_time, 0);
    spikes_comb = real_spikes(data_set, electrode, start_time, end_time, 1);  % combined
    
    % Reject spikes with 0 width
    thres = 10^(-6);

    % Apply threshold to reject reconstructed spikes
    keep_ind = find(c >= thres);
    keep_t = t(keep_ind);
    keep_c = c(keep_ind);
    rej_ind = find(c < thres);
    rej_t = t(rej_ind);
    rej_c = c(rej_ind);

    % Compare spikes
    [tp, fp, fn, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(keep_t, spikes_comb, tolerance);

    scale = 1/(get_dt()*10);

    % Plot results
    fig = figure;
    hold on;
    [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
    %dat_orig = -filter_generic(dat, args);  % Highpass only
    dat = -hp_handle(dat, args);  % Highpass and lowpass
    %p0 = plot(time, dat_orig, 'c');
    p1 = plot(time, dat, 'c');  % Real data
    %stem(t, c*scale, 'k');
    p2 = stem(keep_t, keep_c*scale, 'k');  % Reconstructed spikes
    %stem(rej_t, rej_c*scale, 'Color', [0.5, 0.5, 0.5]);  % Rejected spikes
    p3 = plot(t1, ones(size(t1))*-25, 'b.');  % Sampling times
    p4 = plot(spikes{1}, ones(size(spikes{1}))*-15, 'k.');  % Real spikes
    if length(spikes) > 1
        plot(spikes{2}, ones(size(spikes{2}))*-15, 'kx');  % Real spikes
    end
    %plot([time(1), time(end)], scale*[thres, thres], 'g--');  % Threshold

    % Plot false positives, etc.
    p5 = plot([-1 keep_t(tp)], [0 ; ones(size(tp))*-55], 'g.');
    p6 = plot([-1 keep_t(fp)], [0 ; ones(size(fp))*-45], 'm.');
    p7 = plot([-1 ; spikes_comb(fn)], [0 ; ones(size(fn))*-35], 'r.');

    xlim([time(1), time(end)]);

    leg = legend([p1 p2 p3 p4 p5 p6 p7], 'data', 'reconstructed spikes', 'sampling times', 'real spikes', 'true positives', 'false positives', 'false negatives');
    set(leg, 'Location', 'SouthOutside');
    xlabel('time (sec)');

    % Output error rates
    title({['method: ' method_name ', data set ' int2str(data_set) ', electrode ' int2str(electrode)], ['true positive rate: ' num2str(tp_rate)], ['false positive rate: ' num2str(fp_rate)], ['miss rate: ' num2str(fn_rate)], ['error rate: ' num2str(error_rate)]});

end