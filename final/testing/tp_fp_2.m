function [fn_vec fp_vec] = tp_fp_2(method_name, paramID, data_set, electrode)
% Outputs data for plot of missed spikes (false negative rate) vs false
% positive rate

% Parameters
cond = 'raw';
duration = 10;
tolerance = 0.0025;

[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = data_start;
end_time = start_time + duration;

% get uni params
[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);

% Reconstruct spikes
[t1, result] = run_generic(data_set, [electrode], [start_time], [end_time], hp_handle, diode_handle, method_name, args, start_time, T, window);
[t, c, t1, sigma] = extract_result_single(t1, result);
t = t{1};
c = c{1};
t1 = [t1 t1(end)+T];
sigma = sigma{1};

% 2 methods for getting real spikes: 'real' and threshold

% Real spikes
spikes = real_spikes(data_set, electrode, start_time, end_time, 1);

mult_vec = -2 : 0.25 : 2;
fn_vec = zeros(size(mult_vec));
fp_vec = zeros(size(mult_vec));
for i = 1 : length(mult_vec)

    mult = mult_vec(i);
    
    % Choose threshold based on multiple of std dev  TODO
    thres = mean(c) + mult * std(c);
    %thres = mult * std(c);
    %thres = median(c);

    % Apply threshold to reject reconstructed spikes
    keep_ind = find(c >= thres);
    keep_t = t(keep_ind);
    keep_c = c(keep_ind);
    rej_ind = find(c < thres);
    rej_t = t(rej_ind);
    rej_c = c(rej_ind);

    % Compare spikes
    [tp, fp, fn, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(keep_t, spikes, tolerance);
    
    fn_vec(i) = fn_rate*100;
    fp_vec(i) = length(fp)/duration;%fp_rate;
    
    % PLOT
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
    plot([-1 keep_t(tp)], [0 ; ones(size(tp))*-55], 'g.');
    plot([-1 keep_t(fp)], [0 ; ones(size(fp))*-45], 'm.');
    plot([-1 ; spikes(fn)], [0 ; ones(size(fn))*-35], 'r.');

    xlim([time(1), time(end)]);

    leg = legend('data', 'reconstructed spikes', 'rejected spikes', 'sampling times', 'real spikes', 'threshold', 'sigma', 'samples', 'model samples', 'true positives', 'false positives', 'false negatives');
    set(leg, 'Location', 'SouthOutside');
    xlabel('time (sec)');

    % Output error rates
    title({['method: ' method_name ', data set ' int2str(data_set) ', electrode ' int2str(electrode)], ['true positive rate: ' num2str(tp_rate)], ['false positive rate: ' num2str(fp_rate)], ['miss rate: ' num2str(fn_rate)], ['error rate: ' num2str(error_rate)]});
    %}
end