function [fp_vec, fn_vec] = optthres_reconstruct(data_set, electrodes, method_name, args, hp_handle, diode_handle, T, window)
% Optimizes threshold using training data and returns results (fp/fn rates) on real data

[training_start, training_end, data_start, data_end] = data_division(data_set);

% Time interval for real run
start_time = training_start;  % Run on training data
end_time = start_time + 5;

% ih filter
%[A, B] = butter_hp(400, 10);
%[A, B] = butter_bp(300, 6000, 2);
%args.A_filter = A;
%args.B_filter = B;
cond = 'raw';
%hp_handle = @filtfilt_generic;
%diode_handle = @ideal_diode;

%{
% succ int
T = 0.02;
args.L = 3;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'int';
window = 0.005;
%}

%{
% sinc
T = 0.02;
args.B = 300;
args.L = 10;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'ker';
args.h = @normalized_sinc;
window = 0.1;
%}

%{
% RC
T = 0.01;
args.RC = 0.003;
args.L = 10;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'ker';
args.h = @RC_lowpass;
window = 0.1;
%}

tolerance = 0.01;

fp_vec = [];
fn_vec = [];

% Estimate optimal thresholds
est_rej_thres_gen(data_set, electrodes, cond, hp_handle, diode_handle, method_name, T, window, args, tolerance);

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

    % Get optimal threshold
    thres = choose_threshold(data_set, electrode, method_name, 1);  % Must load, compute using est_rej_thres.m

    % Apply threshold to reject reconstructed spikes
    keep_ind = find(c >= thres);
    keep_t = t(keep_ind);
    keep_c = c(keep_ind);
    rej_ind = find(c < thres);
    rej_t = t(rej_ind);
    rej_c = c(rej_ind);

    % Compare spikes
    [tp, fp, fn, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(keep_t, spikes, tolerance);

    fp_vec = [fp_vec, fp_rate];
    fn_vec = [fn_vec, fn_rate];

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

    %%{

    % Plot false positives, etc.
    plot(keep_t(tp), ones(size(tp))*-55, 'g.');
    plot(keep_t(fp), ones(size(fp))*-45, 'm.');
    plot(spikes(fn), ones(size(fn))*-35, 'r.');

    leg = legend('data', 'reconstructed spikes', 'rejected spikes', 'sampling times', 'real spikes', 'threshold', 'sigma', 'true positives', 'false positives', 'false negatives');
    set(leg, 'Location', 'SouthOutside');
    xlabel('time (sec)');

    % Output error rates
    title({['data set ' int2str(data_set) ', electrode ' int2str(electrode)], ['true positive rate: ' num2str(tp_rate)], ['false positive rate: ' num2str(fp_rate)], ['miss rate: ' num2str(fn_rate)], ['error rate: ' num2str(error_rate)]});

    %}

    %saveas(fig, ['~/Desktop/spk_cmp_' int2str(data_set) '-' int2str(electrode) '.pdf']);

    %%}

end
