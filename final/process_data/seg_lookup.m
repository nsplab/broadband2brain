function num = seg_lookup(data_set, time)
% Gives the segment (file) number of the data segment containing time

filename = [root_dir() 'process_data/data_lookup.mat'];  % Must be the same as in config_lookup.m
load(filename);  % Load table, start_time, end_time

if time < start_time(data_set)
    disp('WARNING: time out of bounds (less than start_time)');
    num = 1;
else if time > end_time(data_set)
    disp('WANRING: time out of bounds (greater than end_time)');
    num = length(table{data_set});
else
    ind = find(table{data_set} <= time);  % Indices of segment files whose start times are <= time
    num = ind(end);  % Index of last segment file whose start time is <= time
end

end
