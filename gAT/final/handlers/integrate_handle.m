function [res next_data] = integrate_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Returns integral of the interval - no diode applied

% Inputs:
% seg_start - segment start time
% seg_end - segment end time
% t1 - interval start time
% T - length of interval
% time - time vector
% dat - data vector
% args = []

% Outputs:
% res - res(1, :) = time, res(2, :) = data

res = zeros(1, 1);
dat_int = integrate(dat, get_dt());
res(1, 1) = dat_int(end);

next_data = 0;

end
