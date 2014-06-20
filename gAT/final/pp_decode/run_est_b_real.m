function [b_0, b_1, b_2] = run_est_b_real(data_set, electrode_num)
% Runs est_b.m using real spikes so that results are saved in .mat file

% Parameters
[d filename] = data_file(data_set, 0, 0);
%filename = 'matthew20080721';
%start_time = 100;
%end_time = 2000;
[start_time, end_time, unused_1, unused_2] = data_division(data_set);  % Training data

elecs = channels_to_use(data_set);
elec = elecs(electrode_num);

% Load data
load(filename);
time_ind = find(time >= start_time - 1 & time < end_time + 1);
time = time(time_ind);
vx = vx(time_ind);
vy = vy(time_ind);

%spike_times = spk(neuron_num).data;
spike_times = real_spikes(data_set, elec, start_time, end_time, 1);  % Combined
spike_ind = find(spike_times >= start_time - 1 & spike_times < end_time + 1);
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

%bs = zeros(length(delay_vec), 3);

delay_vec = get_delay_real(data_set);
delay = delay_vec(electrode_num);

%for l = 1 : length(delay_vec),

    %round(delay_vec(l) / dt)
    spike_v = spike_vector(index - round(delay / dt));

    [b, dev] = glmfit([vx vy],spike_v,'poisson','const','on');

    %bs(l, :) = b;
    %dev_vec(l) = dev;

%end

b_0_new = b(1) - log(get_dt_info());
b_1_new = b(2);
b_2_new = b(3);

b_filename = [root_dir() 'pp_decode/b_temp.mat'];  % TODO: TEMP
load(b_filename);

b_0(electrode_num) = b_0_new;
b_1(electrode_num) = b_1_new;
b_2(electrode_num) = b_2_new;

save(b_filename, 'b_0', 'b_1', 'b_2', '-append');

end
