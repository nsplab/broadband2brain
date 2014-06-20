function [] = config_lookup(data_sets)
% Configures easy lookup of data segment given start time
% Only needs to be run once
% Saves results in a .mat file
% data_sets is a vector of data_set numbers, i.e. [1, 2]

if nargin < 1
    data_sets = [1, 2, 4:19];
end

filename = [root_dir() 'process_data/data_lookup.mat'];

table = {};

% Start/end time for entire data set
start_time = [];
end_time = [];

for i = data_sets
    table{i} = zeros(num_data_segs(i), 1);
    for j = 1 : length(table{i})
        load(data_file(i, 0, j));
        table{i}(j) = time(1);
    end
    start_time(i) = table{i}(1);
    end_time(i) = time(end);
end

save(filename, 'table', 'start_time', 'end_time');
