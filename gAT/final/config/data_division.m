function [training_start, training_end, data_start, data_end] = data_division(data_set)
% Returns start/end of training/data segments
% Values returned are in seconds

% Parameters
offset = 3;  % Buffer in at beginning and end (seconds)
sep = 3;  % Buffer in between (seconds)

[start_time, end_time] = get_total_time(data_set);

training_start = start_time + offset;
data_end = end_time - offset;
training_end = (start_time + end_time - sep) / 2;
data_start = (start_time + end_time + sep) / 2;

end
