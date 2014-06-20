function [x_new, time_new] = cut_segment(time, x, start_time, end_time)
% Prunes x to interval [start_time, end_time]

ind = find(time >= start_time & time < end_time);
x_new = x(ind);
time_new = time(ind);

end
