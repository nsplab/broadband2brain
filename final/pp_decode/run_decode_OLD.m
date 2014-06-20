function [av_RMSE_real av_RMSE_reconstructed] = run_decode(paramID)
% Runs decodes
% Loads spike times
% Applies rejection threshold
% Loads delay, b params

% WARNING: old version

%close all;
%clear all;

% Params
method_name = 'ker';
%paramID = 1;
data_set = 1;
elecs = channels_to_use(data_set);
num_reaches = 200;
plot_flag = 1;
mult = 2;

% TESTING: remove channel 3
%elecs(elecs == 3) = [];

% Time segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start;
end_time = data_end;

% Get real spikes
for i = 1 : length(elecs)
    spike_times(1, i).data = real_spikes(data_set, elecs(i), start_time, end_time, 1);
end

% Get reconstructed spikes, apply rejection threshold
filename = [root_dir() 'reconstruct/saved_spikes/spk_' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
load(filename);  % Loads t, c

for i = 1 : length(elecs)

    elec = elecs(i);

    % Prune reconstructed spikes
    t_e = t{elec};
    c_e = c{elec};
    ind = find(t_e >= start_time & t_e < end_time);
    t_e = t_e(ind);
    c_e = c_e(ind);

    % Apply threshold
    thres = mult * std(c_e);
    ind = find(c_e >= thres);
    t_keep = t_e(ind);
    c_keep = c_e(ind);

    spike_times(2, i).data = t_keep;

end

% Get delay, b
b_file = [root_dir() 'estimate/delay_b/' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
load(b_file);  % Loads params

for i = 1 : length(elecs)
    elec = elecs(i);
    s = params{elec};
    delay_vec(1, i) = s.delay_real;
    delay_vec(2, i) = s.delay;
    b(1, i, :) = s.b_real;
    b(2, i, :) = s.b;
end

% Run decodes
method_names = {'real'; method_name};
[RMSE_px, RMSE_py, RMSE_vx, RMSE_vy] = pp_filter_rand(data_set, start_time, end_time, num_reaches, spike_times, plot_flag, method_names, delay_vec, b);

av_RMSE_real = mean(RMSE_px(1, :)) + mean(RMSE_py(1, :))
av_RMSE_reconstructed = mean(RMSE_px(2, :)) + mean(RMSE_py(2, :))

diff_x = RMSE_px(2, :) - RMSE_px(1, :);
diff_y = RMSE_py(2, :) - RMSE_py(1, :);

% Stats
frac_x = length(find(diff_x > 0)) / length(diff_x);
frac_y = length(find(diff_y > 0)) / length(diff_y);
alpha = 0.05;
[h, p, ci_x] = ttest(RMSE_px(2, :), RMSE_px(1, :), alpha);
[h, p, ci_y] = ttest(RMSE_py(2, :), RMSE_py(1, :), alpha);

% TESTING: close reach plots
close all;

% Plot distribution
figure;
subplot(2, 1, 1);
hold on;
low = -10;
high = 10;
bins = 20;
[x, y freq] = make_histogram(low, high, bins, diff_x);
plot(x, y, 'b');
plot([ci_x(1) ci_x(1)], ylim, 'g--');
plot([ci_x(2) ci_x(2)], ylim, 'g--');
legend('RMSE x', 'confidence interval');
title({['data set: ' int2str(data_set)], 'distribution of (RMSE reconstructed) - (RMSE real)'});
subplot(2, 1, 2);
hold on;
[x, y freq] = make_histogram(low, high, bins, diff_y);
plot(x, y, 'r');
plot([ci_y(1) ci_y(1)], ylim, 'g--');
plot([ci_y(2) ci_y(2)], ylim, 'g--');
legend('RMSE y', 'confidence interval');
disp(['fraction of positive (bad) differences:  x: ' num2str(frac_x) ', y: ' num2str(frac_y)]);
