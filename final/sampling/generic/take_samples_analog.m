function [y sample_times] = take_samples_analog(time, dat, t1, T, args)
% Pass data through comparator that outputs 0 and 1 based on threshold
% Then take successive integrals

try
    thres = args.thres;
catch
    thres = get_thres(args.data_set, args.channel);
end
    
sample_times = [0 T/2];  % whatever

%{
figure;
hold on;
plot(time, dat, 'c');

% Lowpass filter (info spreading)
[A, B] = butter_bp(10, 1000, 2);
args.A_filter = A;
args.B_filter = B;
dat = filter_generic(dat, args);  % Use causal version

plot(time, dat, 'r');
%}

% Prune
ind = find_fast(time, t1, t1+T, get_dt());
dat = dat(ind);
time = time(ind);

% Comparator
dat = dat >= thres;

% Successive integrals
y = succ_int(dat, median(diff(time)), 2);  % slow!

% TESTING: Compute integral of spike
if y(1) > 0
    t = t1 + T - y(2)/y(1);
    w = 0.001;
    ind = find_fast(time, t-w, t+w, get_dt());
    y(3) = get_dt() * sum(dat(ind));
else
    y(3) = 0;
end

% Limit precision
% TODO

%figure;
%hold on;
%plot(time(ind), dat);

end
