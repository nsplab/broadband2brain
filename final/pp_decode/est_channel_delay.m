function [delay_vec, dev_vec, best_delay] = delay_est(data_set, electrode_index)
% Estimates delay between spike times and arm movement

% Author: Alex Wein, July 2011

close all;

% Parameters
[d filename] = data_file(data_set, 0, 0);
%filename = 'matthew20080721';
%start_time = 100;
%end_time = 2000;
[start_time, end_time, unused_1, unused_2] = data_division(data_set);  % Training data

delay_vec = -0.1 : 0.01 : 0.3;  % -200 to 200 ms in intervals of 25 ms
dev_vec = zeros(1, length(delay_vec));

elecs = channels_to_use(data_set);
elec = elecs(electrode_index);

% Load data
load(filename);
time_ind = find(time >= start_time - 1.1*max(delay_vec) & time < end_time + 1.1*max(delay_vec));
time = time(time_ind);
vx = vx(time_ind);
vy = vy(time_ind);

%spike_times = spk(neuron_num).data;
spike_times = real_spikes(data_set, elec, start_time, end_time, 1);  % Combined
spike_ind = find(spike_times >= start_time - 1.1*max(delay_vec) & spike_times < end_time + 1.1*max(delay_vec));
spike_times = spike_times(spike_ind);

dt = (max(time) - min(time)) / (length(time) - 1);

% Generate spike_vector from spike_times
spike_vector = zeros(1, length(time));
i = 1;  % Iterates over time indices
j = 1;  % Iterates over spike time indices

while j <= length(spike_times) && spike_times(j) < time(1),
    j = j + 1;
end
while i < length(time) && j <= length(spike_times),
    if spike_times(j) >= time(i) && spike_times(j) < time(i+1),
        spike_vector(i) = 1;
        while j <= length(spike_times) && spike_times(j) < time(i+1),
            j = j + 1;
        end
    end
    i = i + 1;
end

index = find(time >= start_time & time < end_time);
vx = vx(index);
vy = vy(index);

bs = zeros(length(delay_vec), 3);

for l = 1 : length(delay_vec),

    disp([int2str(l) '/' int2str(length(delay_vec))]);

    %round(delay_vec(l) / dt)
    spike_v = spike_vector(index - round(delay_vec(l) / dt));

    [b, dev] = glmfit([vx vy],spike_v,'poisson','const','on');

    bs(l, :) = b;
    dev_vec(l) = dev;

end

[best_dev best_ind] = min(dev_vec);
best_delay = delay_vec(best_ind(1));
bs(best_ind(1), :)
disp(['optimal delay: ' num2str(best_delay) ' sec']);

figure;
plot(delay_vec, dev_vec);
title({['deviance vs delay (electrode ' num2str(elec) ')'],['optimal delay = ' num2str(best_delay) ' sec']});
xlabel('delay (seconds)');
ylabel('deviance');
