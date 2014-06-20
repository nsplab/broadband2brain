% Figure 3 - "ROC" curves and error rate vs T
% This script plots multiple ROC curves (different channels) for some fixed
% T, tolerance, delta (refractory period)

close all;
clear all;

% Parameters
data_set = 1;
channels = channels_to_use(data_set);
channels_to_plot{1} = [16 15 13];%channels;%
channels_to_plot{2} = [3 15 12];%channels;%
channels_to_plot{3} = [3 15 12];%channels;%
duration = 60;
start_offset = 60;
delta = 0.0011;  % refractory period
tolerance = 0.005;  % for deciding if a spike is "correct"
T = 0.100;

% Data segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start + start_offset;
end_time = start_time + duration;

% Plot
fig = figure;
hold on;

methods = {'Analog Thresholding', 'gAT', 'twodelta'};
colors = {[1 0 0], [0 0 1], [0 1 0]};

vals{1} = [];
vals{2} = [];
vals{3} = [];

for i = 1 : length(channels)
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

    num_thres = 50;
    thres_vec = linspace(1*mean(dat)+0*max(dat), 0.2*mean(dat)+0.8*max(dat), num_thres);  % thresholds to try

    for m = 1 : length(methods)

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

                    t = t_TD(c ~= 0); % Need to remove times corresponding to invalid spikes
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
        
        %methods{m}
        %electrode
        %val2
        disp([methods{m} ' channel ' int2str(electrode) ', error rate = ' num2str(val2)]);
        vals{m} = [vals{m} val2];
        
        if any(channels_to_plot{m} == electrode)
            [fn_roc fp_roc] = convexify(fn_roc, fp_roc);  % convexify
            plot_errorbar(fn_roc, fp_roc, fp_roc, fp_roc, colors{m}, 1);
            plot_dots(fn_roc(ind2), fp_roc(ind2), darken(colors{m}));
        end
        %title(['ROC: T = ' num2str(T)]);
        %title(['ROC: T = ' num2str(T) ', best thres = ' num2str(thres_vec(ind)) ', error rate = ' num2str(val/length(real_spk))]);
        p1 = plot_errorbar(-10, -10, -10, -10, colors{1}, 1);
        p2 = plot_errorbar(-10, -10, -10, -10, colors{2}, 1);
        p3 = plot_errorbar(-10, -10, -10, -10, colors{3}, 1);
        leg = legend([p1 p2 p3], 'AT', 'gAT', 'twodelta');
        set(leg, 'FontSize', 28);
        xlim([0 1]);
        ylim([0 1]);
        
    end
    
end

sort(vals{1})
sort(vals{2})
sort(vals{3})

xlabel('False Negatives per True Spike', 'FontSize', 28);
ylabel('False Positives per True Spike', 'FontSize', 28);
set(gca, 'FontSize', 24);
name = 'median';  % all/median
dir = [root_dir() '../paper2plots/final/fig3/'];
saveas(fig, [dir 'ROC_' name '_' num2str(1000*T)], 'epsc');
saveas(fig, [dir 'ROC_' name '_' num2str(1000*T)], 'fig');
