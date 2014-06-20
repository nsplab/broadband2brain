function [t] = run_threshold(data_set, electrodes, start_time, end_time, cond)
% Runs analog thresholding

% Parameters
window = 0;
T = 50;  % Length of itervals, can't exceed memory

[t1, result] = iterate_intervals(data_set, electrodes, [start_time], [end_time], cond, start_time, T, window, @run_threshold_handle, []);
t1 = t1{1};

% Extract t from results
t = cell(length(electrodes), 1);
for i = 1 : length(result)
    t{i} = result{i}(1, :);
end

end
