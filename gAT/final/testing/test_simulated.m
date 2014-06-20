% Tests generic reconstruction method on simulated data
% This allows accurate comparison of real/reconstructed amplitudes
% Mostly copy/pasted from run_generic, run_generic_handle
% Bypasses iterate_itervals

close all;
clear all;

% Parameters
method_name = 'RMSE'; %'ker';
paramID = 14;
duration = 1;  % in seconds (same dt as real data)

% get uni params
[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);
% hp_handle and diode_handle go unused - we assume simulated data already
% conditioned

% TESTING: modify params
args.numiter = 10;
%args.delta_t = get_dt()*10000;

% make sure duration is a multiple of T
duration = floor(duration / T) * T;

lambda = 15;  % spikes per second

% generate simulated spike train
time = -2*window : get_dt() : duration + 2*window;  % some unnecessary buffer
prob = lambda * get_dt();  % probability of a spike at each vector entry
dat = unifrnd(0, 1, 1, length(time)) < prob;
ind = find(time < 0 | time >= duration);
dat(ind) = zeros(size(ind));  % no spikes in window

t_real = time(find(dat ~= 0));  % real spike times
c_real = normrnd(0.03, 0.006, 1, length(t_real));
dat = +dat;  % convert logical (0,1) vector to numeric vector
dat(find(dat ~= 0)) = c_real / get_dt();

% add noise
delta_t_ind = 5;  % real delta_t in units of dt
delta_t_real = delta_t_ind * get_dt();
sigma_real = 0.0001 / get_dt(); %0.01;
ind = find(dat == 0);
noise = max(0, normrnd(0, sigma_real, 1, length(ind)));
i = 1;
while(i <= length(noise))  % step function noise
    j = 0;
    while(j < delta_t_ind && i + j <= length(noise))
        noise(i + j) = noise(i);
        j = j + 1;
    end
    i = i + delta_t_ind;
end
dat(ind) = noise;

% Setup function
setup_func = eval(['@setup_' method_name]);
setup_obj = setup_func(T, window, args);

% Append to args
args.setup_obj = setup_obj;
args.hp_handle = hp_handle;
args.diode_handle = diode_handle;
args.method_name = method_name;

% aggregated vectors
y_tot = [];
time_y_tot = [];
t_tot = [];
c_tot = [];
y_model = [];

next_data = struct();

sigma_sum = 0;

% reconstruct, interval (of size T) by interval
for t1 = 0 : T : duration - T/2  % T/2 instead of T in case numerical issues
    
    % get data
    ind = find(time >= t1 - window & time < t1 + T + window);
    
    % take samples
    sample_func = eval(['@take_samples_' args.method_name]);
    [y, sample_times] = sample_func(time, dat(ind), t1, T, args);
    
    y_tot = [y_tot y'];
    time_y_tot = [time_y_tot sample_times];
    
    % run algorithm
    rec_func = eval(['@reconstruct_' args.method_name]);
    [t, c, sigma, elapsed_time] = rec_func(y, T, args, next_data);
    
    sigma_sum = sigma_sum + sigma;
    
    % model-reconstructed samples
    [y_m] = model_samples(t1, sample_times, T, args, next_data, t+t1, c, sigma);
    y_model = [y_model y_m'];
    
    % Pass spikes on to next interval
    next_data = struct();
    next_data.t = t;
    next_data.c = c;

    % Limit spike times to 0, T]
    ind = find(t >= 0 & t < T);
    t = t(ind, 1);
    c = c(ind, 1);
    
    t = t + t1;
    
    t_tot = [t_tot t'];
    c_tot = [c_tot c'];
    
end

t = t_tot;
c = c_tot;

% stats
sigma_real
sigma_av = sigma_sum / length(0 : T : duration - T/2)
c_real_av = mean(c_real)
c_av = mean(c)  % wanring: not fair comparison - no rejection threshold so tiny spikes included
y_real_av = mean(y_tot)
y_model_av = mean(y_model)
delta_t = args.delta_t
delta_t_real

% TODO: threshold

% plot results
figure;
hold on;
plot(time, dat * get_dt(), 'c');
stem(t_real, c_real, 'c');
stem(t, c, 'k');
y_sc = max(c) / (max(y_tot) * 2);
plot(time_y_tot, y_tot * y_sc, 'm.-');
plot(time_y_tot, y_model * y_sc, 'b.-');
t1_vec = 0 : T : duration - T/2;
plot(t1_vec, ones(size(t1_vec))*-0.002, 'b.');  % Computation intervals

xlabel('time (sec)');
title('RMSE-based reconstruction');
legend('data', 'real spikes', 'reconstructed spikes', 'real samples', 'model-reconstructed samples', 'computation intervals');

% plot kernel
%{
sinc_t = 0 : 0.00001 : duration;
sinc_x = zeros(size(sinc_t));
for i = 1 : length(sinc_x)
    sinc_x(i) = args.h(sinc_t(i) - mean(sinc_t), args);
end
plot(sinc_t, sinc_x, 'r');

kernel_integral = sum(sinc_x) * 0.00001
%}

%{

% theoretical samples
y_the = zeros(1, length(time_y_tot));
for l = 1 : length(y_the)
    for k = 1 : length(c_real)
        y_the(l) = y_the(l) + c_real(k) * args.h(time_y_tot(l) - t_real(k), args);
        if l == floor(length(y_the)*4/5)
            alone = c_real(k)
            contr = c_real(k) * args.h(time_y_tot(l) - t_real(k), args)
        end
    end
end
plot(time_y_tot, y_the, 'g.-');

% theoretical samples (2)
y_the = zeros(1, length(time_y_tot));
for l = 1 : length(y_the)
    for i = 1 : length(time)
        y_the(l) = y_the(l) + dat(i) * args.h(time_y_tot(l) - time(i), args) * get_dt();
        if l == floor(length(y_the)*4/5)
            contr = dat(i) * args.h(time_y_tot(l) - time(i), args) * get_dt();
            if contr ~= 0
                alone = dat(i)
                with_dt = dat(i) * get_dt()
                contr
            end
        end
    end
end
plot(time_y_tot, y_the, 'k.-');

%}