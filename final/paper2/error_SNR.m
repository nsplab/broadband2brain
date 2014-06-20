% scatter error rate (false pos + false neg) with SNR
% error rate calculated by finding best point on ROC
% some fixed T value

close all;

% Parameters
data_set = 1;
duration = 60;
T = 0.01;
num_thres = 20;
tolerance = 0.0025;
delta = 0.0011;  % refractory period

colors = {[1 0 0], [0 0 1]};

% Data segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start;
end_time = start_time + duration;
t_block = start_time : T : end_time;
t_block = t_block(1 : end-1);  % t_block contains start times

electrodes = channels_to_use(data_set);
err_vec_AT = zeros(size(electrodes));
err_vec_aFRI = zeros(size(electrodes));

SNR_vec = SNR_test(data_set, electrodes);

for i = 1 : length(electrodes)
    electrode = electrodes(i);

    % Get data
    buff = 0.1;
    [time, dat] = get_data(data_set, electrode, start_time-buff, end_time+buff, 'raw');

    % Highpass
    [A, B] = butter_bp(300, 6000, 2);  % Change filter order?
    args.A_filter = A;
    args.B_filter = B;
    hp_handle = @filter_generic;  % Use causal version
    dat = -hp_handle(dat, args);

    % Get real spikes
    real_spk = real_spikes(data_set, electrode, start_time, end_time, 1);

    thres_vec = linspace(0.9*mean(dat)+0.1*max(dat), 0.3*mean(dat)+0.7*max(dat), num_thres);  % thresholds to try

    fp_roc_AT = zeros(1, length(thres_vec));
    fn_roc_AT = zeros(1, length(thres_vec));
    fp_roc_aFRI = zeros(1, length(thres_vec));
    fn_roc_aFRI = zeros(1, length(thres_vec));

    % Loop through thres
    for j = 1 : length(thres_vec)
        thres = thres_vec(j);

        t_AT = [];
        t_aFRI = [];

        % Loop through intervals
        for k = 1 : length(t_block)

            % get data segment
            seg_ind = find_fast(time, t_block(k), t_block(k)+T, get_dt());
            seg = dat(seg_ind);

            % get reconstructed spikes in interval
            % analog thresholding
            num_rec_spikes_AT = max(seg) > thres;
            if num_rec_spikes_AT == 1 && (length(t_AT) < 1 || T/2+t_block(k) > t_AT(end)+delta)
                t_AT = [t_AT T/2+t_block(k)];
            end
            % aFRI
            seg_comp = seg > thres;
            y = succ_int(seg_comp, get_dt(), 2);
            [t_F, c] = reconstruct_analog(y, T, [], []);
            if c > 0 && (length(t_aFRI) < 1 || t_F+t_block(k) > t_aFRI(end)+delta)
                t_aFRI = [t_aFRI t_F+t_block(k)];
            end

        end
        
        [~, ~, ~, ~, ~, fp_rate, fn_rate] = compare_spikes(t_AT, real_spk, tolerance);
        fp_roc_AT(j) = fp_rate;
        fn_roc_AT(j) = fn_rate;

        [~, ~, ~, ~, ~, fp_rate, fn_rate] = compare_spikes(t_aFRI, real_spk, tolerance);
        fp_roc_aFRI(j) = fp_rate;
        fn_roc_aFRI(j) = fn_rate;
        
    end

    % For this channel, find best threshold
    [val ind_AT] = min(fp_roc_AT + fn_roc_AT);
    err_vec_AT(i) = val;
    [val ind_aFRI] = min(fp_roc_aFRI + fn_roc_AT);
    err_vec_aFRI(i) = val;
    
    % Plot ROC curve using raw compare_spikes
    figure;
    hold on;
    plot_errorbar(fn_roc_AT, fp_roc_AT, fp_roc_AT, fp_roc_AT, colors{1}, 1);
    plot_errorbar(fn_roc_aFRI, fp_roc_aFRI, fp_roc_aFRI, fp_roc_aFRI, colors{2}, 1);
    plot_dots(fn_roc_AT(ind_AT), fp_roc_AT(ind_AT), [0 1 0]);
    plot_dots(fn_roc_aFRI(ind_aFRI), fp_roc_aFRI(ind_aFRI), [0 1 0]);
    p1 = plot_errorbar(-10, -10, -10, -10, colors{1}, 1);
    p2 = plot_errorbar(-10, -10, -10, -10, colors{2}, 1);
    legend([p1 p2], 'AT', 'aFRI');
    title(['channel ' int2str(electrode)]);
    xlim([0 1]);
    ylim([0 1]);
    
end

figure;
hold on;
scatter(SNR_vec, err_vec_AT, 'r');
scatter(SNR_vec, err_vec_aFRI, 'b');
legend('AT', 'aFRI');
xlabel('SNR');
ylabel('error rate (false pos + false neg per real spike)');
title(['T = ' num2str(T)]);