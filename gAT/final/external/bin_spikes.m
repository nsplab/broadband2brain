function spike_vector = bin_spikes(time, spike_times)
% Converts spike times to vector of spike counts
% spike_times(i) is number of spikes in interval [ time(i) , time(i+1) )

% Inputs:
% time - time vector
% spike_times - list of spike times (in order)

spike_vector = zeros(1, length(time)-1);
i = 1;  % Iterates over time indices
j = 1;  % Iterates over spike time indices

while j <= length(spike_times) && spike_times(j) < time(1),
    j = j + 1;
end
while i < length(time) && j <= length(spike_times),
    if spike_times(j) >= time(i) && spike_times(j) < time(i+1),
        %spike_vector(i) = spike_vector(i) + 1;
        while j <= length(spike_times) && spike_times(j) < time(i+1),
            j = j + 1;
            spike_vector(i) = spike_vector(i) + 1;
        end
    end
    i = i + 1;
end
