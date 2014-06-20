function [y sample_times] = take_samples_AT(time, dat, t1, T, args)
% Returns a single sample, 0 or 1 (for whether threshold was crossed)

try
    thres = args.thres;
catch
    thres = get_thres(args.data_set, args.channel);
end
    
sample_times = [0];  % whatever

% Prune
ind = find_fast(time, t1, t1+T, get_dt());
dat = dat(ind);
%time = time(ind);

% Comparator
y = max(dat) >= thres;

end
