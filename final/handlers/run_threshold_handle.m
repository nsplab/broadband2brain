function [res next_data] = run_threshold_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Analog thresholding
% res(1, :) = times

% Get eyeballed threshold from config/get_threshold.m
threshold = get_threshold(data_set, elec, 'ih');

% Run reconstruction
t = analog_threshold(time, dat, threshold);

% Limit spike times to [seg_start, seg_end]
ind = find(t >= seg_start & t < seg_end);
t = t(ind, 1);

% Output
res = zeros(1, length(t));

res(1, :) = t';

next_data = 0;

end
