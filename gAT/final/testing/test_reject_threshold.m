% Plots missed spikes and false positives verses threshold used to ignore reconstructed spikes
% Uses training data

close all;

% Parameters
data_set = 1;
real_load_b = 0;
FRI_load_b = 1;
real_b_type = 'real';
FRI_b_type = 'real';
cond = 'ih';
t0 = 0;
T = 0.02;
L = 5;
K = 1;
delta_t = get_dt();
numiter = 3;
t_resolution = T / 100;

elec = 16;
tolerance = 0.005;

[training_start, training_end, data_start, data_end] = data_division(data_set);

% TESTING
training_end = training_start + 20;

% Get real spikes
t_real = real_spikes(data_set, elec, training_start, training_end, 1);

% Reconstruct spikes using FRI_int_diode (in data and also in training if needed)
[t, c, t1] = run_FRI_int_diode(data_set, [elec], training_start, training_end, cond, t0, T, L, K, delta_t, numiter, t_resolution);
t = t{1};
c = c{1};

thres_to_try = 0 : 0.001 : 0.06;
fpos = zeros(size(thres_to_try));
fneg = zeros(size(thres_to_try));

for i = 1 : length(thres_to_try)
    thres = thres_to_try(i);
    ind = find(c >= thres);
    [tp, fp, fn] = compare_spikes(t(ind), t_real, tolerance);
    fpos(i) = length(fp);
    fneg(i) = length(fn);
end

% Plot
figure;
hold on;
plot(thres_to_try, fpos, 'b');
plot(thres_to_try, fneg, 'r');
plot(thres_to_try, fpos + fneg, 'k');
legend('false positives', 'missed spikes', 'total errors');
xlabel('threshold');

[val, in] = min(fpos + fneg);
disp(['best threshold = ' num2str(thres_to_try(in))]);
