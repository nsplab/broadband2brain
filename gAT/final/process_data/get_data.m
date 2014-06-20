function [time, dat] = get_data(data_set, electrodes, start_time, end_time, cond)
% Returns time and data vectors

% Inputs:
% data_set
% electrodes - vector of electrode numbers
% start_time in seconds
% end_time in seconds
% cond - data conditioning type

% Outputs:
% time vector (row)
% dat(elec_num, :) = data vector

[t1, result] = iterate_intervals(data_set, electrodes, [start_time], [end_time], cond, start_time, end_time-start_time, 0, @get_data_handle, []);
time = zeros(1, size(result{1, 1}, 2));
dat = zeros(length(electrodes), size(result{1, 1}, 2));

for elec_num = 1 : length(electrodes)
    if elec_num == 1
        time(1, :) = result{elec_num, 1}(1, :);
    end
    dat(elec_num, :) = result{elec_num, 1}(2, :);
end

end
