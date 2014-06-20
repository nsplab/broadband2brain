% ACTUALLY FIG 3

% Figure 2 - "ROC" curves and error rate vs T
% Plot error rate vs T for some fixed
% tolerance, delta (refractory period)

close all;
clear all;

% Parameters
data_set = 2;
channels = channels_to_use(data_set);
%channels_to_plot{1} = [16 15 13];
%channels_to_plot{2} = [3 15 12];
duration = 60;
start_offset = 60;
delta = 0.0011;  % refractory period
tolerance = 0.005;  % for deciding if a spike is "correct"
T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods

debug = 0;  % if 1 plot individual ROC curves
load = 0;

if load == 0

    % Data segment
    [training_start, training_end, data_start, data_end] = data_division(data_set);
    start_time = data_start + start_offset;
    end_time = start_time + duration;

    % Plot
    if debug
        figure;
        hold on;
    end

    methods = {'Analog Thresholding', 'gAT', 'twodelta'};
    colors = {[1 0 0], [0 0 1], [0 1 0]};

    data_AT = cell(1, 16);  % data_AT{channel} = vector of error_rates (parallel to T_vec)
    data_aFRI = cell(1, 16);
    data_TD = cell(1, 16);
    thres_AT = cell(1, 16);
    thres_aFRI = cell(1, 16);
    thres_TD = cell(1, 16);

    for i = 1 : length(channels)
        fprintf(sprintf('i: %d / %d\n', i, length(channels)));
        electrode = channels(i);

        % Get data
        buff = 0.1;
        [time, dat] = get_data(data_set, electrode, start_time-buff, end_time+buff, 'raw');

        %figure;
        %plot(time, dat);
        %return;

        % Highpass
        [A, B] = butter_bp(300, 6000, 2);  % Change filter order?
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        dat = -hp_handle(dat, args);

        % Get real spikes
        real_spk = real_spikes(data_set, electrode, start_time, end_time, 1);
        %av_spk_rate = length(real_spk)/duration

        num_thres = 30;
        thres_vec = linspace(1*mean(dat)+0*max(dat), 0.2*mean(dat)+0.8*max(dat), num_thres);  % thresholds to try

        for m = 1 : length(methods)
            fprintf(sprintf('m: %d / %d\n', m, length(methods)));

            for p = 1 : length(T_vec)
                fprintf(sprintf('p: %d / %d\n', p, length(T_vec)));
                T = T_vec(p);

                t_block = start_time : T : end_time;
                t_block = t_block(1 : end-1);  % t_block contains start times

                % Loop through thres
                for j = 1 : length(thres_vec)
                    thres = thres_vec(j);

                    real_ind = 1;  % index in real_spk

                    t_rec = [];

                    % Loop through intervals
                    for k = 1 : length(t_block)

                        % find real spikes in interval
                        num_real_spikes = 0;
                        t_real = [];
                        while real_ind <= length(real_spk) && real_spk(real_ind) <= t_block(k)+T
                            t_real = [t_real real_spk(real_ind) - t_block(k)];
                            num_real_spikes = num_real_spikes + 1;
                            real_ind = real_ind + 1;
                        end

                        % get data segment
                        seg_ind = find_fast(time, t_block(k), t_block(k)+T, get_dt());
                        seg = dat(seg_ind);

                        % get reconstructed spikes in interval
                        if m == 1
                            % analog thresholding
                            t_AT = T/2;
                            num_rec_spikes_AT = max(seg) > thres;

                            t = t_AT;
                            num_rec_spikes = num_rec_spikes_AT;
                        elseif m == 2
                            % aFRI
                            seg_comp = seg > thres;
                            y = succ_int(seg_comp, get_dt(), 2);
                            [t_F, c] = reconstruct_analog(y, T, [], []);
                            num_rec_spikes_F = (c ~= 0);

                            t = t_F;
                            num_rec_spikes = num_rec_spikes_F;
                        elseif m == 3
                            % twodelta
                            seg_comp = seg > thres;
                            y = succ_int(seg_comp, get_dt(), 4);
                            [t_TD, c] = reconstruct_twodelta(y, T, [], []);
                            num_rec_spikes_TD = sum(c ~= 0);

                            t = t_TD;
                            num_rec_spikes = num_rec_spikes_TD;
                        end

                        % Check that there are actually reconstructed spikes, then
                        %   check that this is either the first one or that the refractory period has passed
                        if num_rec_spikes > 0 && (length(t_rec) < 1 || t_block(k)+min(t) > t_rec(end)+delta)
                            t_rec = [t_rec t_block(k)+t'];
                        end

                    end

                    % for ROC
                    [~, ~, ~, ~, ~, fp_rate, fn_rate] = compare_spikes(t_rec, real_spk, tolerance);
                    fp_roc(j) = fp_rate;
                    fn_roc(j) = fn_rate;

                end

                % Plot ROC curve using raw compare_spikes
                [val2 ind2] = min(fp_roc + fn_roc);
                best_thres = thres_vec(ind2);
                if m == 1
                    data_AT{electrode}(p) = val2;
                    thres_AT{electrode}(p) = best_thres;
                elseif m == 2
                    data_aFRI{electrode}(p) = val2;
                    thres_aFRI{electrode}(p) = best_thres;
                elseif m == 3
                    data_TD{electrode}(p) = val2;
                    thres_TD{electrode}(p) = best_thres;
                end

                if debug
                    if nnz(channels_to_plot{m} == electrode) > 0
                        [fn_roc fp_roc] = convexify(fn_roc, fp_roc);  % convexify
                        plot_errorbar(fn_roc, fp_roc, fp_roc, fp_roc, colors{m}, 1);
                        plot_dots(fn_roc(ind2), fp_roc(ind2), darken(colors{m}));
                    end
                    title(['ROC: T = ' num2str(T)]);
                    %title(['ROC: T = ' num2str(T) ', best thres = ' num2str(thres_vec(ind)) ', error rate = ' num2str(val/length(real_spk))]);
                    p1 = plot_errorbar(-10, -10, -10, -10, colors{1}, 1);
                    p2 = plot_errorbar(-10, -10, -10, -10, colors{2}, 1);
                    p3 = plot_errorbar(-10, -10, -10, -10, colors{3}, 1);
                    legend([p1 p2 p3], 'AT', 'gAT', 'twodelta');
                    xlim([0 1]);
                    ylim([0 1]);
                end
            end

        end

    end

    save([root_dir() 'paper2/mat/fig3b' int2str(data_set)]);

end





    
clearvars -except data_set;
load([root_dir() 'paper2/mat/fig3b' int2str(data_set)]);

dir = [root_dir() '../paper2plots/final/fig3/'];

% All channels (separate curves)
fig = figure;
hold on;
for electrode = channels
    p1 = plot_errorbar(1000*T_vec, data_AT{electrode}, 0, 0, [1 0 0], 1);
    p2 = plot_errorbar(1000*T_vec, data_aFRI{electrode}, 0, 0, [0 0 1], 1);
    p3 = plot_errorbar(1000*T_vec, data_TD{electrode}, 0, 0, [0 1 0], 1);
end
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Errors per True Spike', 'FontSize', 28);
%title('All Channels');
leg = legend([p1 p2 p3], 'AT', 'gAT', 'twodelta');
set(leg, 'Location', 'NorthWest', 'FontSize', 28);
ylim([0 1.2]);
xlim(1000*[0 0.052]);
set(gca, 'FontSize', 24);
saveas(fig, [dir 'err_all' int2str(data_set)], 'epsc');
saveas(fig, [dir 'err_all' int2str(data_set)], 'fig');

% Bootci errorbars
fig = figure;
hold on;
mean_AT = zeros(1, length(T_vec));
ci_low_AT = zeros(1, length(T_vec));
ci_hi_AT = zeros(1, length(T_vec));
mean_aFRI = zeros(1, length(T_vec));
ci_low_aFRI = zeros(1, length(T_vec));
ci_hi_aFRI = zeros(1, length(T_vec));
mean_TD = zeros(1, length(T_vec));
ci_low_TD = zeros(1, length(T_vec));
ci_hi_TD = zeros(1, length(T_vec));
for i = 1 : length(T_vec)
    vec_AT = [];
    vec_aFRI = [];
    vec_TD = [];
    for elec = channels
        vec_AT = [vec_AT data_AT{elec}(i)];
        vec_aFRI = [vec_aFRI data_aFRI{elec}(i)];
        vec_TD = [vec_TD data_TD{elec}(i)];
    end
    mean_AT(i) = mean(vec_AT);
    mean_aFRI(i) = mean(vec_aFRI);
    mean_TD(i) = mean(vec_TD);
    ci = bootci(1000, @mean, vec_AT);
    ci_low_AT(i) = ci(1);
    ci_hi_AT(i) = ci(2);
    ci = bootci(1000, @mean, vec_aFRI);
    ci_low_aFRI(i) = ci(1);
    ci_hi_aFRI(i) = ci(2);
    ci = bootci(1000, @mean, vec_aFRI);
    ci_low_TD(i) = ci(1);
    ci_hi_TD(i) = ci(2);
end
p1 = plot_errorbar(1000*T_vec, mean_AT, ci_low_AT, ci_hi_AT, [1 0 0]);
p2 = plot_errorbar(1000*T_vec, mean_aFRI, ci_low_aFRI, ci_hi_aFRI, [0 0 1]);
p3 = plot_errorbar(1000*T_vec, mean_TD, ci_low_TD, ci_hi_TD, [0 1 0]);
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Errors per True Spike', 'FontSize', 28);
%title('Mean with bootci confidence intervals for the mean');
leg = legend([p1 p2 p3], 'AT', 'gAT', 'twodelta');
set(leg, 'Location', 'NorthWest', 'FontSize', 28);
ylim([0 1.2]);
xlim(1000*[0 0.052]);
set(gca, 'FontSize', 24);
saveas(fig, [dir 'err_mean' int2str(data_set)], 'epsc');
saveas(fig, [dir 'err_mean' int2str(data_set)], 'fig');

% Thresholds
for electrode = channels
    figure;
    hold on;
    p1 = plot_errorbar(T_vec, thres_AT{electrode}, 0, 0, [1 0 0], 1);
    p2 = plot_errorbar(T_vec, thres_aFRI{electrode}, 0, 0, [0 0 1], 1);
    p3 = plot_errorbar(T_vec, thres_aFRI{electrode}, 0, 0, [0 1 0], 1);
    xlabel('Sampling Period (T)');
    ylabel('Best Threshold');
    leg = legend([p1 p2 p3], 'AT', 'gAT', 'twodelta');
    set(leg, 'Location', 'NorthWest');
    title(['Channel ' int2str(electrode) ', best thres (aFRI, T = 0.005) = ' num2str(thres_aFRI{electrode}(3))]);
    fprintf(['Electrode #' int2str(electrode) ':\n']);
    for i = 1:length(T_vec)
        fprintf(['    T = ' num2str(T_vec(i)) ':\t' num2str(thres_AT{electrode}(i), 6) '\t' num2str(thres_aFRI{electrode}(i), 6) '\t' num2str(thres_TD{electrode}(i), 6) '\n'])
    end
end

