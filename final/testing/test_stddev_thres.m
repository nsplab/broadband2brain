function [sig, best_multiple] = test_stddev_thres(method_name, paramID, data_set, elec)
% Rejection threshold chosen by multiple of std dev
% Plots error rates vs multiple

%close all;
%clear all;

start_time = 100;
end_time = 200;

if nargin == 0
    method_name = 'int';
    paramID = 1;
    data_set = 1;
    elec = 9;
end

tolerance = 0.01;

% Load reconstructed spikes
filename = [root_dir() 'reconstruct/saved_spikes/spk_' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
load(filename);  % Loads t, c

% Load real spikes
t_real = real_spikes(data_set, elec, start_time, end_time, 1);

t = t{elec};
c = c{elec};
ind = find(t >= start_time & t < end_time);
t = t(ind);
c = c(ind);

sig = std(c);

As = 0 : 0.1 : 5;

fp = zeros(size(As));
fn = zeros(size(As));

for i = 1 : length(As)

    ind = find(c >= sig * As(i));
    [true_positives, false_positives, false_negatives, error_rate, tp_rate, fp_rate, fn_rate] = compare_spikes(t(ind), t_real, tolerance);
    fp(i) = fp_rate;
    fn(i) = fn_rate;

end

figure;
hold on;
plot(As, fp, 'b');
plot(As, fn, 'r');
plot(As, fp + fn, 'k');
legend('false positive rate', 'false negative rate', 'error rate');
xlabel('multiple');
title(['method: ' method_name ', data set: ' int2str(data_set) ', elec: ' int2str(elec)]);

[val, best_ind] = min(fp + fn);
std_dev = sig
best_multiple = As(best_ind)
