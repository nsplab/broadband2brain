function [time, data_vec] = load_data_file(data_set, electrodes, file_num, cond)
% Loads a data file

% Inputs:
% data_set
% electrodes - list of electrodes
% file_num - which data segment to load
% cond - if data should be conditioned (see config/data_file.m)

% Outputs:
% time - time vector (row)
% data_vec - data_vec(elec, :) = data for electrode elec

for i = 1 : length(electrodes)
    elec = electrodes(i);
    load(data_file(data_set, elec, file_num, cond));  % Loads time and dat
    if i == 1  % First iteration
        data_vec = zeros(length(electrodes), length(dat));
    end
    data_vec(i, :) = dat;
end

time = time';  % Return row vector

end
