function spike_vector = build_spike_vec(time, spike_times)
% Converts spike times to vector of 0's and 1's (or higher numbers in rare cases)
% spike_vector(i) = 1 iff there is a spike between time(i) and time(i+1)
% Processes a single list of spike times

% Inputs:
% time - time vector
% spike_times - list of spike times

spike_vector = zeros(1, length(time));
i = 1;  % Iterates over time indices
j = 1;  % Iterates over spike time indices

while j <= length(spike_times) && spike_times(j) < time(1),
    j = j + 1;
end
while i < length(time) && j <= length(spike_times),
    if spike_times(j) >= time(i) && spike_times(j) < time(i+1),
        while j <= length(spike_times) && spike_times(j) < time(i+1),
            j = j + 1;
            spike_vector(i) = spike_vector(i) + 1;
        end
    end
    i = i + 1;
end
