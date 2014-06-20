function spike_vector = spike_times_to_vec(time, spike_times)
% Converts spike times to vector of 0's and 1's
% Processes all neurons

% Inputs:
% time - time vector
% spike_times - has form spike_times(neuron_num).data = vector of spike times

C = length(spike_times);
spike_vector = zeros(C, length(time));
for c = 1 : C,
    i = 1;  % Iterates over time indices
    j = 1;  % Iterates over spike time indices

    while j <= length(spike_times(c).data) && spike_times(c).data(j) < time(1),
        j = j + 1;
    end
    while i < length(time) && j <= length(spike_times(c).data),
        if spike_times(c).data(j) >= time(i) && spike_times(c).data(j) < time(i+1),
            while j <= length(spike_times(c).data) && spike_times(c).data(j) < time(i+1),
                spike_vector(c, i) = spike_vector(c, i) + 1;
                j = j + 1;
            end
        end
        i = i + 1;
    end
end
