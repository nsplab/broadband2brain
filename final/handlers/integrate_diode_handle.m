function [res next_data] = integrate_diode_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Same as integrate_handle but applies ideal diode

res = integrate_handle(data_set, elec, seg_start, seg_end, t1, T, time, ideal_diode(dat), args);

next_data = 0;

end
