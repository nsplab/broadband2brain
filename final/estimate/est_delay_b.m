function [real_dir_vec dir_vec dir_err dir_err_vec] = est_delay_b(method_name, paramIDs, data_sets, period, electrodes)
% Estimates delay and b params given training data spike times
% Plots results
% For decodes, need to run this on training data = first half of
% non-force-field data
% For tuning curves (nemo channel 5) need to run this on all data,
% separately for baseline, force field, washout

method_name

tic

%close all;
%clear all;

% Parameters
%method_name = 'int';
%paramID = 1;
%data_set = 1;
%elec = 16;
mult = 0;  % Multiple of sigma for choosing rejection threshold

%methods = {'int', 'ker', 'RMSE'};

if nargin < 3
    data_sets = [1, 2];  % nemo channel 5
end

% TESTING: pref dir plot
real_dir_vec = [];
dir_vec = [];
real_dep_vec = [];
dep_vec = [];

%for method_ind = 2; %[1, 2];
%    method_name = methods{method_ind};
    for paramID = paramIDs
        for data_set = data_sets;
            data_set
            if nargin < 5
                electrodes = channels_to_use(data_set);
            end
            for elec = electrodes; %channels_to_use(data_set)

                min_delay = -0.1;  % Must be >= 0
                max_delay = 0.9;
                delta_delay = 0.05;

                if nargin < 4 || period < 2
                    % Training data
                    [training_start, training_end, data_start, data_end] = data_division(data_set);
                    start_time = training_start;
                    end_time = training_end;
                else
                    % Baseline/ff/washout (depending on period argument)
                    [start_time, end_time, ff_start, ff_end, wash_start, wash_end] = get_total_time(data_set);
                    if period == 2
                        start_time = ff_start;
                        end_time = ff_end;
                    elseif period == 3
                        start_time = wash_start;
                        end_time = wash_end;
                    end
                end

                % Load reconstructed spikes
                filename = [root_dir() 'reconstruct/saved_spikes/spk_' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
                load(filename);  % Loads t, c

                % Prune reconstructed spikes
                t = t{elec};
                c = c{elec};
                %original_t_starttime = t(1)
                ind = find(t >= start_time & t < end_time);
                t = t(ind);
                c = c(ind);
                
                % Load real spikes
                t_real = real_spikes(data_set, elec, start_time, end_time, 1);
                t_real_orig = t_real;  % for plot
                
                start_time
                rec_spike_start = t(1)
                real_spike_start = t_real_orig(1)
                
                if strcmp(method_name, 'gibbs')
    
                    % Choose threshold based on multiple of std dev  TODO
                    %thres = mean(c) + mult * std(c);
                    %thres = mult * std(c);
                    %thres = median(c);
                    
                    % Alternative thresholding: median-based
                    c_sorted = sort(c);
                    frac = 3/4;
                    thres = c_sorted(round(length(c) * frac));
                    
                elseif strcmp(method_name, 'analog') || strcmp(method_name, 'AT') || strcmp(method_name, 'twodelta')
        
                    thres = 0.00001;

                else

                    % Choose threshold for 30 spikes per second
                    num_to_keep = round(30 * (end_time-start_time));
                    c_sorted = sort(c);
                    thres = c_sorted(length(c_sorted) - num_to_keep);

                end
                
                % Reject spikes below thres
                ind = find(c >= thres);
                t_keep = t(ind);
                c_keep = c(ind);
                t_keep_orig = t_keep;  % for plot
                c_keep_orig = c_keep;
                
                % Load arm velocity
                [d filename] = data_file(data_set, 0, 0);
                load(filename);  % Loads time, vx, vy
                
                %size_time = size(time)
                %min_time = min(time)
                %max_time = max(time)

                %%%

                % only use reaches as training data
                %%{
                mothr = 5; % movement onset threshold (cm/s)
                wnd = [0.1 0.5]; % window relative to movement onset (s)
                %load([monkey session '.mat'],'time','cio','cof','tio','tof','vx','vy');
                fs = round(1/median(diff(time))); % sampling rate of prepared data (Hz)
                swnd = round(wnd*fs);
                indcus = find(diff(tof)>=1); % indices of 'cue'
                indgos = find(cio==1 & [0;diff(cof)==-1]); % indices of 'go'
                indres = find(tio==1 & [0;diff(tof)<=-1]); % indices of 'rew'
                tmove = zeros(length(indres),2); % trials x [move start, move end] in seconds
                for itrial = 1:length(indres) % for each rewarded trial
                    indre = indres(itrial);
                    indgo = indgos(find(indgos < indre)); indgo = indgo(end);
                    speed = sqrt(vx(indgo:indre).^2+vy(indgo:indre).^2);
                    spthr = find(speed > mothr); indmo = indgo + spthr(1) - 1; % movement onset
                    tmove(itrial,:) = time(indmo+swnd);
                end

                % Prune tmove to start_time, end_time
                ind = find(tmove(:, 1) >= start_time);
                tmove = tmove(ind, :);
                ind = find(tmove(:, 2) < end_time);
                tmove = tmove(ind, :);

                % Stitch together reaches
                time_new = [];
                vx_new = [];
                vy_new = [];
                last_time = time(1);
                t_real_new = [];
                t_keep_new = [];
                c_keep_new = [];
                t_keep = t_keep';
                c_keep = c_keep';

                for i = 1 : size(tmove, 1)
                    sec_start = tmove(i, 1);
                    sec_end = tmove(i, 2);
                    
                    if sec_start > sec_end
                        disp('WARNING: sec_start > sec_end (est_delay_b.m)');
                        continue;
                    end
                    
                    ind = find(time >= sec_start & time < sec_end);
                    change_t = last_time - time(ind(1)) + get_dt_info();
                    time_new = [time_new ; time(ind) + change_t];
                    last_time = time_new(end);
                    vx_new = [vx_new ; vx(ind)];
                    vy_new = [vy_new ; vy(ind)];

                    t_real_new = [t_real_new ; t_real(find(t_real >= sec_start & t_real < sec_end)) + change_t];
                    ind = find(t_keep >= sec_start & t_keep < sec_end);
                    t_keep_new = [t_keep_new ; t_keep(ind) + change_t];
                    c_keep_new = [c_keep_new ; c_keep(ind)];
                end

                time = time_new;
                vx = vx_new;
                vy = vy_new;
                t_real = t_real_new;
                t_keep = t_keep_new;
                c_keep = c_keep_new;

                %%}

                % OLD
                %{
                % Prune arm velocity
                ind = find(time >= start_time & time < end_time);
                time = time(ind);
                vx = vx(ind);
                vy = vy(ind);
                %}

                %%%


                % Build spike vectors
                spk_vec = build_spike_vec(time, t_keep);
                spk_vec_real = build_spike_vec(time, t_real);

                delay_step_ind = round(delta_delay / get_dt_info());
                min_delay_ind = round(min_delay / (delay_step_ind * get_dt_info())) * delay_step_ind;
                max_delay_ind = round(max_delay / (delay_step_ind * get_dt_info())) * delay_step_ind;

                spk_ind_adj = 1 - min_delay_ind : length(spk_vec) - max_delay_ind;
                
                % Vectors to plot
                shift_vec = min_delay_ind : delay_step_ind : max_delay_ind;
                delay_vec = shift_vec * get_dt_info();
                dev_vec = zeros(size(delay_vec));
                dev_vec_real = zeros(size(delay_vec));
                b_vec = zeros(length(delay_vec), 3);
                b_vec_real = zeros(size(b_vec));
                b_vec_er = zeros(size(b_vec));
                b_vec_real_er = zeros(size(b_vec));

                for i = 1 : length(shift_vec)

                    disp([int2str(i) ' of ' int2str(length(shift_vec))]);

                    shift = shift_vec(i);

                    time_ind_adj = spk_ind_adj(1) + shift : spk_ind_adj(end) + shift;
                    
                    % glmfit for reconstructed spikes
                    [b, dev, stats] = glmfit([vx(time_ind_adj) vy(time_ind_adj)], spk_vec(spk_ind_adj), 'poisson', 'const', 'on');
                    b(1) = b(1) - log(get_dt_info());  % Correct for the fact that glmfit's lambda is spikes per timestep (not per second)
                    dev_vec(i) = dev;
                    b_vec(i, :) = b;
                    b_vec_er(i, :) = stats.se;

                    % glmfit for real spikes
                    [b, dev, stats] = glmfit([vx(time_ind_adj) vy(time_ind_adj)], spk_vec_real(spk_ind_adj), 'poisson', 'const', 'on');
                    b(1) = b(1) - log(get_dt_info());  % Correct for the fact that glmfit's lambda is spikes per timestep (not per second)
                    dev_vec_real(i) = dev;
                    b_vec_real(i, :) = b;
                    b_vec_real_er(i, :) = stats.se;
                    
                    % TESTING: plot spikes to check for offset problem
                    % Everything seems fine...
                    %{
                    figure;
                    stem(time_ind_adj, spk_vec(spk_ind_adj));
                    title('test 1');
                    figure;
                    stem(time_ind_adj, spk_vec_real(spk_ind_adj));
                    title('test 2');
                    %}

                end

                % Get best params
                [val, ind_real] = min(dev_vec_real);
                delay_real = delay_vec(ind_real);
                b_real = b_vec_real(ind_real, :);
                [val, ind] = min(dev_vec);
                delay = delay_vec(ind);
                b = b_vec(ind, :);

                % Get tuning curves
                best_shift_real = shift_vec(ind_real);
                time_ind_adj = spk_ind_adj(1) + best_shift_real : spk_ind_adj(end) + best_shift_real;  % Copied from above
                [theta_real, freq_real, freq_real_b, pref_dir_real] = tuning_curve(vx(time_ind_adj), vy(time_ind_adj), spk_vec_real(spk_ind_adj), b_real);
                best_shift = shift_vec(ind);
                time_ind_adj = spk_ind_adj(1) + best_shift : spk_ind_adj(end) + best_shift;  % Copied from above
                [theta, freq, freq_b, pref_dir] = tuning_curve(vx(time_ind_adj), vy(time_ind_adj), spk_vec(spk_ind_adj), b);

                % For pref dir plots
                real_dir_vec = [real_dir_vec pref_dir_real];
                dir_vec = [dir_vec pref_dir];
                real_dep_vec = [real_dep_vec range(freq_real_b)];
                dep_vec = [dep_vec range(freq_b)];

                % Plot results
                fig = figure;

                % Deviance
                subplot(5, 2, 1);
                plot(delay_vec, dev_vec_real, 'k');
                legend('real', 'Location', 'NorthEast');
                xlabel('delay (sec)');
                title({['method: ' method_name ', data set: ' int2str(data_set) ', elec: ' int2str(elec)], 'deviance'});

                subplot(5, 2, 2);
                plot(delay_vec, dev_vec, 'b');
                legend(method_name, 'Location', 'NorthEast');
                xlabel('delay (sec)');
                title('deviance');

                % b params
                colors = {'b', 'r', 'g'};
                for i = 1 : 3
                    subplot(5, 2, i + 2);
                    hold on;
                    errorbar(delay_vec, b_vec_real(:, i), b_vec_real_er(:, i), 'k');
                    errorbar(delay_vec, b_vec(:, i), b_vec_er(:, i), 'b');
                    %legend('real', 'reconstructed', 'Location', 'Best');
                    xlabel('delay (sec)');
                    title(['b' int2str(i - 1)]);
                end

                %{
                % Spike train
                subplot(3, 2, 6);
                hold on;
                stem(t, c, 'Color', [0.5, 0.5, 0.5]);
                stem(t_keep, c_keep, 'b');
                plot([t(1) t(end)], [thres thres], 'g--');
                stem(t_real, -thres*ones(size(t_real)), 'k');
                %}

                % Tuning curves

                subplot(5, 2, 6);  % Overlay
                hold on;
                plot(theta_real, freq_real / get_dt_info(), 'k');  % Normalize lambda correctly            
                plot(theta_real, freq_real_b, 'k--');
                plot(theta, freq / get_dt_info(), 'b');
                plot(theta, freq_b, 'b--');
                xlabel('theta (radians)');
                ylabel('lambda (spikes/sec)');
                title('tuning curves');
                axis tight;

                subplot(5, 2, 7);  % Real
                hold on;
                plot(theta_real, freq_real / get_dt_info(), 'k');  % Normalize lambda correctly            
                plot(theta_real, freq_real_b, 'k--');
                %plot([pref_dir_real pref_dir_real], ylim, 'k--');
                %plot([pref_dir pref_dir], ylim, 'b--');
                xlabel('theta (radians)');
                ylabel('lambda (spikes/sec)');
                title(['tuning curve: depth ' num2str(range(freq_real_b))]);
                axis tight;

                subplot(5, 2, 8);  % Reconstructed
                hold on;
                plot(theta, freq / get_dt_info(), 'b');
                plot(theta, freq_b, 'b--');
                xlabel('theta (radians)');
                ylabel('lambda (spikes/sec)');
                title(['tuning curve: depth ' num2str(range(freq_b))]);
                axis tight;
                
                % TESTING
                %test = setdiff(t_keep, intersect(t, t_keep))
                
                start_time
                rec_spike_start = t(1)
                real_spike_start = t_real_orig(1)
                
                % Spike train
                subplot(5, 2, [9 10]);
                hold on;
                start_time = t_real_orig(2);
                end_time = start_time + 0.5;
                scale = 1/(get_dt()*10);
                [args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);
                [time, dat] = get_data(data_set, [elec], start_time - window, end_time, 'raw');
                dat = -hp_handle(dat, args);  % Highpass (and flip)
                [dat, time] = cut_segment(time, dat, start_time, end_time);  % Cut out window
                plot(time, dat, 'c');  % Real data
                xlim([start_time, end_time]);
                %stem(t, c*scale, 'k');
                stem(t, c*scale, 'Color', [0.5, 0.5, 0.5]);  % Rejected spikes
                stem(t_keep_orig, -c_keep_orig*scale, 'b');  % Reconstructed spikes
                t1 = start_time - mod(start_time, T) : T : end_time;  % Reverse-engineer t1
                plot(t1, ones(size(t1))*-25, 'r.');  % Sampling times
                plot(t_real_orig, ones(size(t_real_orig))*-15, 'k.');  % Real spikes
                plot([time(1), time(end)], scale*[thres, thres], 'g--');  % Threshold
                %stem(t, sigma, 'm');  % Plot sigma
                xlabel('time (sec)');

                % Save plot
                mkdir([root_dir() 'tuning_curves']);
                saveas(fig, [root_dir() 'tuning_curves/' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID) '.pdf']);

                % Save params in .mat
                clear params;
                b_file = [root_dir() 'estimate/delay_b/' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
                try            
                    load(b_file);
                end
                s.delay_real = delay_real;
                s.b_real = b_real;
                s.delay = delay;
                s.b = b;
                s.thres = thres;
                params{elec} = s;
                save(b_file, 'params');
                
                %close(fig);  % TESTING: close figure
                
                % 2 tuning curves
                %{
                fig = figure;  % Overlay
                hold on;
                %plot_line(theta_real, freq_real / get_dt_info(), 'k');  % Normalize lambda correctly            
                p1 = plot_line(180/pi*theta_real, freq_real_b, [0 0 0]);
                [val ind] = max(freq_real_b);
                plot_dots(180/pi*pref_dir_real, freq_real_b(ind), [0 0 0]);
                %plot_line(theta, freq / get_dt_info(), 'b');
                p2 = plot_line(180/pi*theta, freq_b, [0 0 1]);
                [val ind] = max(freq_b);
                plot_dots(180/pi*pref_dir, freq_b(ind), [0 0 1]);
                xlabel('Angle (degrees)', 'FontSize', 28);
                ylabel('Av. Firing Rate (spikes/sec)', 'FontSize', 28);
                %title('tuning curve reconstruction', 'FontSize', 18);
                axis tight;
                leg = legend([p1 p2], 'True', 'Rec');
                set(leg, 'FontSize', 28, 'Location', 'NorthWest');
                set(gca, 'FontSize', 24);
                
                dir = [root_dir() '../paper2plots/final/fig5/'];
                saveas(fig, [dir 'tuning_' int2str(elec)], 'epsc');
                %}
                
                % save to file instead, used by fig5.m
                dir2 = [root_dir() '/paper2/mat/tuning'];
                save([dir2 '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)], 'theta_real', 'freq_real_b', 'pref_dir_real', 'theta', 'freq_b', 'pref_dir', 'freq_real', 'freq');

            end
        end
    end
%end

% TESTING: pref dir plots
%%{
figure;
hold on;
low = -pi; %min(min(real_dir_vec), min(dir_vec));
high = pi; %max(max(real_dir_vec), max(dir_vec));
plot_line([low, high], [low, high], [0 0 0], 1);
plot_dots(real_dir_vec, dir_vec, [0 0 1]);
title('preferred direction', 'FontSize', 18);
xlabel('real', 'FontSize', 18);
ylabel('gibbs', 'FontSize', 18);
xlim([low high]);
ylim([low high]);

% modulation depth plot
figure;
hold on;
low = 0;
high = 10;
plot_line([low, high], [low, high], [0 0 0], 1);
plot_dots(real_dep_vec, dep_vec, [0 0 1]);
title('depth of modulation', 'FontSize', 18);
xlabel('real', 'FontSize', 18);
ylabel('gibbs', 'FontSize', 18);
xlim([low high]);
ylim([low high]);
%%}

% dir_err: average absolute error in preferred direction
dif = abs(real_dir_vec - dir_vec)
dif = min(dif, 2*pi-dif)
dir_err = mean(dif)

dir_err_vec = dif;

toc
