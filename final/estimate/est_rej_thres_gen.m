function [] = est_rej_thres_gen(data_set, elecs, cond, hp_handle, diode_handle, method_name, T, window, args, tolerance);
% Estimates reject threshold for generic reconstruction and saves it

L = args.L;
K = args.K;
numiter = args.numiter;
t_resolution = args.t_resolution;
id = method_name;

% Parameters
offset = 10;
duration = 60;  % Length of training data to use (seconds)
thres_res = 0.0001;  % Resolution of threshold values to try

[training_start, training_end, data_start, data_end] = data_division(data_set);

start_time = training_start + offset;
end_time = start_time + duration;

% 2 methods of getting real spikes  TODO: choose 1

% Get real spikes
for i = 1 : length(elecs)
    t_real{i} = real_spikes(data_set, elecs(i), start_time, end_time, 1);
end

% Run thresholding
%t_real = run_threshold(data_set, elecs, start_time, end_time, cond);

% Reconstruct spikes
[t1, result] = run_generic(data_set, elecs, start_time, end_time, hp_handle, diode_handle, method_name, args, 0, T, window);
[t, c, t1] = extract_result_single(t1, result);

% Choose threshold
for i = 1 : length(elecs)
    thresholds(i) = choose_threshold(data_set, elecs(i), id, 0, t{i}, c{i}, t_real{i}, tolerance, thres_res);  % Saves thresholds in .mat file
end

end
