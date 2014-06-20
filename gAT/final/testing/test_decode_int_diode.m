% Runs some sample decodes and plots results
% Compares 'real' spikes to FRI_int_diode

close all;

% Parameters
data_set = 1;
real_load_b = 1;
FRI_load_b = 1;
real_b_type = 'real';
FRI_b_type = 'real';
cond = 'ih';
t0 = 0;
T = 0.02;
L = 5;
K = 1;

delta_t = get_dt();
numiter = 3;
t_resolution = T / 100;

load_thres = 0;
thres_id = 'test';
num_reaches = 5;
thres_training_duration = 10;  % Length of training data for choosing thresholds for rejecting FRI spikes
tolerance = 0.005;
thres_res = 0.001;  % Threshold optimization resolution

[training_start, training_end, data_start, data_end] = data_division(data_set);

% Get reach start/end times in data
disp('get reaches');
[reach_start, reach_end] = get_reaches(data_set, data_start, data_end);
reach_start = reach_start(1 : num_reaches);
reach_end = reach_end(1 : num_reaches);
disp('done');

% Get real spikes (in data and also in training if needed)
disp('start get spikes');
elecs = channels_to_use(data_set);
for i = 1 : length(elecs)
    real_spikes_data{i} = real_spikes(data_set, elecs(i), reach_start(1), reach_end(end), 1);
    if real_load_b == 0 || load_thres == 0
        if real_load_b == 0
            end_training_get = training_end;
        else
            end_training_get = training_start + thres_training_duration;
        end
        real_spikes_training{i} = real_spikes(data_set, elecs(i), training_start, end_training_get, 1);
    end
end
disp('done');

% Reconstruct spikes using FRI_int_diode (in data and also in training if needed)
disp('start reconstruct');
[t_data, c_data, t1] = run_FRI_int_diode(data_set, channels_to_use(data_set), reach_start(1), reach_end(end), cond, t0, T, L, K, delta_t, numiter, t_resolution);
if FRI_load_b == 0 || load_thres == 0
    if FRI_load_b == 0
        end_training_get = training_end;
    else
        end_training_get = training_start + thres_training_duration;
    end
    [t_training, c_training, t1] = run_FRI_int_diode(data_set, channels_to_use(data_set), training_start, end_training_get, cond, t0, T, L, K, delta_t, numiter, t_resolution);
end
disp('done');

% Apply threshold to eliminate small spikes (threshold is computed or loaded)
disp('start threshold');
for i = 1 : length(elecs)
    if load_thres == 1
        thres = choose_threshold(data_set, elecs(i), thres_id, load_thres);
    else
        thres = choose_threshold(data_set, elecs(i), thres_id, load_thres, t_training{i}, c_training{i}, real_spikes_training{i}, tolerance, thres_res);
    end
    ind = find(c_data{i} >= thres);
    t_data{i} = t_data{i}(ind);
    c_data{i} = c_data{i}(ind);
end
disp('done');

% Run pp_filter_rand.m
% spike_times(method_num, neuron_num).data = list of spike times
for i = 1 : length(elecs)
    spike_times(1, i).data = real_spikes_data{i};
    spike_times(2, i).data = t_data{i};
end
[RMSE_px, RMSE_py, RMSE_vx, RMSE_vy] = pp_filter_rand(data_set, reach_start(1)-0.01, reach_end(end)+0.01, num_reaches, spike_times, 1, {'analog threshold'; 'successive integral'});
