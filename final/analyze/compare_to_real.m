function [true_pos_rate, false_pos_rate, false_neg_rate] = compare_to_real(data_set, electrodes, hp_handle, diode_handle, method_name, args, T, window)
% Given method and parameters, optimizes rejection threshold and outputs rate of true/false positives/negatives

[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = training_start;  % Run on training data
end_time = start_time + 5;

cond = 'raw';
tolerance = 0.01;

% Estimate optimal thresholds   TODO: only necessary if parameters have changed, can be commented otherwise
%est_rej_thres_gen(data_set, electrodes, cond, hp_handle, diode_handle, method_name, T, window, args, tolerance);
% NOTE: now just choosing best threshold based on this data

true_pos_rate = zeros(size(electrodes));
false_pos_rate = zeros(size(electrodes));
false_neg_rate = zeros(size(electrodes));

for i = 1 : length(electrodes)

    electrode = electrodes(i);

    %ticID = tic;

    % Reconstruct spikes
    [t1, result] = run_generic(data_set, [electrode], [start_time], [end_time], hp_handle, diode_handle, method_name, args, start_time, T, window);
    [t, c, t1, sigma] = extract_result_single(t1, result);
    t = t{1};
    c = c{1};
    t1 = [t1 t1(end)+T];
    sigma = sigma{1};

    %reconstruct_time = toc(ticID)

    % 2 methods for getting real spikes: 'real' and threshold

    %ticID = tic;

    % Real spikes
    spikes = real_spikes(data_set, electrode, start_time, end_time, 1);

    %real_time = toc(ticID)

    % Analog threshold spikes
    %spikes = run_threshold(data_set, [electrode], start_time, end_time, cond);
    %spikes = spikes{1};

    %ticID = tic;

    % Get optimal threshold
    %thres = choose_threshold(data_set, electrode, method_name, 1);  % Must load, compute using est_rej_thres.m
    % Just choose best one for this data
    thres_res = 0.0001;  % Not actually used anymore
    thres = choose_threshold(data_set, electrode, 'opt', 0, t, c, spikes, tolerance, thres_res);

    %choose_thres_time = toc(ticID)
    %ticID = tic;

    % Apply threshold to reject reconstructed spikes
    keep_ind = find(c >= thres);
    keep_t = t(keep_ind);
    keep_c = c(keep_ind);
    rej_ind = find(c < thres);
    rej_t = t(rej_ind);
    rej_c = c(rej_ind);

    %apply_time = toc(ticID)
    %ticID = tic;

    % Compare spikes
    [tp, fp, fn, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(keep_t, spikes, tolerance);

    %compare_time = toc(ticID)

    true_pos_rate(i) = tp_rate;
    false_pos_rate(i) = fp_rate;
    false_neg_rate(i) = fn_rate;


    % TESTING: plot

    %{

    scale = 1/(get_dt()*10);

    % Plot results
    fig = figure;
    hold on;
    [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
    dat = -hp_handle(dat, args);  % Highpass
    plot(time, dat, 'c');  % Real data
    %stem(t, c*scale, 'k');    
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

    %saveas(fig, ['~/Desktop/spk_cmp_' int2str(data_set) '-' int2str(electrode) '.pdf']);

    %}

end
