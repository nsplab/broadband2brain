function t = analog_threshold(time, z, threshold, delta)
% Spike detection using a threshold

% Inputs:
% time - time vector
% z - data
% threshold
% delta - minimum time allowed between spikes

% Outputs:
% t - spike times (column vector)

if nargin < 4
    delta = 0.0008;
end

N = length(z);

t = []; % Reconstructed spike times
j = 1; % Index in t

% Find spikes
lastTime = -delta; % Last spike time
for i = 1 : N,
    if abs(z(i)) >= threshold && time(i) >= lastTime + delta,  % TODO: is analog threshold 2-sided?
        t(j, 1) = time(i);
        lastTime = time(i);
        j = j + 1;
    end
end

end
