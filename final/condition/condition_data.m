function [] = condition_data(data_sets, highpass_func, highpass_args, diode_func, diode_args, cond, window, flip)
% Prepares data for sampling
% Applies flip, highpass filter and/or diode
% Re-saves data

% Inputs:
% data_sets - vector of data set numbers to process
% highpass_func - function handler or 0 for none
% highpass_args - arguments to above
% diode_func - function handler or 0 for none
% diode_args - arguments to above
% cond - folder name for new data files
% window - [before, after] how much surrounding data is needed for highpass to avoid edge effects (seconds)
% flip - 0 or 1 for whether data should be negated (to make spikes positive)

if nargin < 7
    flip = 1;
end
if nargin < 8
    window = [2, 2];
end

if strcmp(cond, 'raw')
    return;
end

% Make folder
mkdir([root_dir 'data'], cond)

% Convert window to units of indices
window = round(window / get_dt());

prev_time = [];
cur_time = [];
next_time = [];
prev_data = [];
cur_data = [];
next_data = [];

for data_set = data_sets

    % Initialize
    el = channels_to_use(data_set);
    num_el = length(el);
    prev_time = [];
    prev_data = [];
    [cur_time, cur_data] = load_data_file(data_set, el, 1, 'raw');
    [next_time, next_data] = load_data_file(data_set, el, 2, 'raw');

    for seg = 1 : num_data_segs(data_set)

        % Get current segment (with windows)
        time = cur_time';
        data = cur_data;
        orig_start = 1;
        orig_end = size(data, 2);
        if size(prev_time, 1) > 0
            data = [prev_data(:, end-window(1)+1 : end) data];
            orig_start = orig_start + window(1);
            orig_end = orig_end + window(1);
        end
        if size(next_time, 1) > 0
            data = [data next_data(:, 1 : window(2))];
        end
        % Condition and save
        if flip
            data = -data;  % Flip
        end
        for i = 1 : num_el
            dat = data(i, :)';
            if ~isnumeric(highpass_func)
                dat = highpass_func(dat, highpass_args);  % Highpass
            end
            if ~isnumeric(diode_func)
                dat = diode_func(dat, diode_args);  % Diode
            end
            dat = dat(orig_start : orig_end, 1);
            save(data_file(data_set, el(i), seg, cond), 'time', 'dat');  % Save
        end

        % Update time/data vectors
        prev_time = cur_time;
        prev_data = cur_data;
        cur_time = next_time;
        cur_data = next_data;
        if seg >= num_data_segs(data_set) - 1
            next_time = [];
            next_data = [];
        else
            [next_time, next_data] = load_data_file(data_set, el, seg+2, 'raw');
        end

    end
end 

end
