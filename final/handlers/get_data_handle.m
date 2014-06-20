function [res next_data] = get_data_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Simple handler that returns time and data for plotting

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

if time(1) < seg_start || time(end) >= seg_end
    ind = find(time >= seg_start & time < seg_end);
    time = time(ind);
    dat = dat(ind);
end

res = zeros(2, length(time));
res(1, :) = time;
res(2, :) = dat;

next_data = 0;

end
