function [samples, samples_model, tp_rate, fp_rate] = test_generic_reconstruct(paramID, mult)
% Runs generic reconstruct method
% Plots spike train
% NOTE: return value assumes a single electrode is being analyzed

%close all;
%clear all;

% Parameters
data_set = 1;
electrodes = [11]; %channels_to_use(data_set);
cond = 'raw';
duration = 60; %%%%%

method_name = 'analog'; %'gibbs' for conference
if nargin < 1
    paramID = 15; %34 is 1 kHz (for gibbs)
end

[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = data_start;
end_time = start_time + duration;

% get uni params instead of above
[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);
if nargin < 2
    mult = 0;
end

% threshold for analog
% if not specified, uses defaults in paper2/get_thres.m
%args.thres = 75.5388;%57.4095;

% TESTING: modify params
%%args.numiter = 50;
%%args.mean_iter = 10;
%%args.t_resolution = args.t_resolution / 10;
%diode_handle = @fullwave_rectify;

%args.delta_t = get_dt() * 1000;
%{
T = 0.03;
args.L = 30;
args.K = 2;
args.numiter = 10;
args.delta_t = 0.0005; %get_dt() * 15;
args.s1 = -1000; %-16000;
args.s2 = -800; %-20000;
%args.t_resolution = args.t_resolution / 10;
% Filter
%%{
[A, B] = butter_hp(1000, 2);
%[A, B] = butter_bp(300, 6000, 2);
args.A_filter = A;
args.B_filter = B;
hp_handle = @filter_generic;  % Causal version
%hp_handle = @filtfilt_generic;  % Noncausal version
%%}
%}

% TESTING: clear file for saved samples
filename = [root_dir() 'mat_files/samples.mat'];
t_samples = [];
samples = [];
samples_model = [];
runtime = [];
save(filename, 't_samples', 'samples', 'samples_model', 'runtime');

tolerance = 0.0025;

% Estimate optimal thresholds   TODO: only necessary if parameters have changed, can be commented otherwise
%est_rej_thres_gen(data_set, electrodes, cond, hp_handle, diode_handle, method_name, T, window, args, tolerance);

for electrode = electrodes

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

    % Analog threshold spikes
    %spikes = run_threshold(data_set, [electrode], start_time, end_time, cond);
    %spikes = spikes{1};

    % Get optimal threshold (old)
    %thres = choose_threshold(data_set, electrode, method_name, 1);  % Must load, compute using est_rej_thres.m (very outdated)
    
    if strcmp(method_name, 'gibbs')
    
        % Choose threshold based on multiple of std dev  TODO
        thres = mean(c) + mult * std(c);
        %thres = mult * std(c);
        %thres = median(c);
        
    elseif strcmp(method_name, 'analog') || strcmp(method_name, 'AT')
        
        thres = 0.00001;
        
    else
    
        % Choose threshold for 30 spikes per second
        num_to_keep = round(30 * duration);
        c_sorted = sort(c);
        thres = c_sorted(length(c_sorted) - num_to_keep);

    end
    
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
    dat = -hp_handle(dat, args);  % Highpass
    plot(time, dat, 'c');  % Real data
    %stem(t, c*scale, 'k');    
    stem(keep_t, keep_c*scale, 'k');  % Reconstructed spikes
    stem(rej_t, rej_c*scale, 'Color', [0.5, 0.5, 0.5]);  % Rejected spikes
    plot(t1, ones(size(t1))*-25, 'b.');  % Sampling times
    plot(spikes, ones(size(spikes))*-15, 'k.');  % Real spikes
    plot([time(1), time(end)], scale*[thres, thres], 'g--');  % Threshold
    stem(t, sigma, 'm');  % Plot sigma
    
    % Samples
    load(filename);
    plot(t_samples, samples, 'm.-');
    plot(t_samples, samples_model, 'b.-');

    %%{

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

    %saveas(fig, ['~/Desktop/spk_train/spk_cmp_' method_name '_' int2str(data_set) '_' int2str(electrode) '_' int2str(paramID) '.pdf']);
    
    %%}
    
    % TESTING: plot for conference proceeding
    fig = figure;
    hold on;
    [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
    dat = -hp_handle(dat, args);  % Highpass
    time_new = 1000*(time-time(1));
    p1 = plot(time_new, dat, 'c');  % Real data
    %stem(t, c*scale, 'k');    
    p2 = stem(1000*(keep_t-time(1)), keep_c*scale/2, 'k');  % Reconstructed spikes
    set(p2, 'LineWidth', 2);
    %stem(rej_t, rej_c*scale, 'Color', [0.5, 0.5, 0.5]);  % Rejected spikes
    %plot(t1, ones(size(t1))*-25, 'b.');  % Sampling times
    h = plot(1000*(spikes-time(1)), ones(size(spikes))*-30, 'kx');  % Real spikes
    set(h, 'LineWidth', 2, 'MarkerSize', 10);
    %plot([time(1), time(end)], scale*[thres, thres], 'g--');  % Threshold
    %stem(t, sigma, 'm');  % Plot sigma
    
    save([root_dir() '../data'], 'time', 'dat');
    
    % Samples
    %load(filename);
    %plot(t_samples, samples*2, 'm.-');
    %plot(t_samples, samples_model*2, 'b.-');

    %%{

    % Plot false positives, etc.
    %plot([-1 keep_t(tp)], [0 ; ones(size(tp))*-55], 'g.');
    %plot([-1 keep_t(fp)], [0 ; ones(size(fp))*-45], 'm.');
    %plot([-1 ; spikes(fn)], [0 ; ones(size(fn))*-35], 'r.');

    xlim([time_new(1), time_new(end)]);
    
    method_title = 'annihilating filter';
    leg = legend([p1 p2 h], 'data', method_title, 'true spikes');
    set(leg, 'Location', 'SouthOutside', 'FontSize', 14);
    xlabel('time (ms)', 'FontSize', 14);
    ylabel('amplitude (mV)', 'FontSize', 14);
    
    % Output error rates
    %title({['method: ' method_name ', data set ' int2str(data_set) ', electrode ' int2str(electrode)], ['true positive rate: ' num2str(tp_rate)], ['false positive rate: ' num2str(fp_rate)], ['miss rate: ' num2str(fn_rate)], ['error rate: ' num2str(error_rate)]}, 'FontSize', 12);



end