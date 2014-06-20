function [y sample_times] = take_samples_twodelta(time, x, t1, T, args)
% Pass data through comparator that outputs 0 and 1 based on threshold
% Then take successive integrals

try
    thres = args.thres;
catch
    thres = get_thres(args.data_set, args.channel);
end

sample_times = 0; % This isn't actually used in the code (to the best of my knowledge)

% Prune
ind = find_fast(time, t1, t1+T, get_dt());
x = x(ind);
time = time(ind);

% Comparator
x = x >= thres;

% Successive integrals
y = succ_int(x, median(diff(time)), 4);  % slow!

end
