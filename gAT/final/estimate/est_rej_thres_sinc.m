function [] = est_rej_thres(data_set, elecs, id, cond, T, L, K, B, delta_t, numiter, t_resolution, tolerance);
% Estimates reject threshold for FRI_sinc_diode and saves it

% Parameters
duration = 60;  % Length of training data to use (seconds)
thres_res = 0.0001;  % Resolution of threshold values to try

[training_start, training_end, data_start, data_end] = data_division(data_set);

% 2 methods of getting real spikes  TODO: choose 1

% Get real spikes
%for i = 1 : length(elecs)
%    t_real{i} = real_spikes(data_set, elecs(i), training_start, training_start + duration, 1);
%end

% Run thresholding
t_real = run_threshold(data_set, elecs, training_start, training_start + duration, cond);

% Get FRI_int_diode spikes
[t, c, t1] = run_FRI_sinc_diode(data_set, elecs, training_start, training_start + duration, cond, 0, T, L, K, B, delta_t, numiter, t_resolution);

% Choose threshold
for i = 1 : length(elecs)
    thresholds(i) = choose_threshold(data_set, elecs(i), id, 0, t{i}, c{i}, t_real{i}, tolerance, thres_res);  % Saves thresholds in .mat file
end

end
