function [x, y freq] = make_histogram(low, high, bins, data)
% Outputs x and y so that plot(x, y) generates a normalized histogram
vec = linspace(low, high, bins+1);
freq = zeros(bins, 1);
for i = 1 : length(freq),
    freq(i) = length(find(data >= vec(i) & data < vec(i+1)));
end
freq = freq / (sum(freq) * (vec(2)-vec(1)));  % Normalize
[x, y] = plot_histogram(vec, freq);
freq = freq * (vec(2) - vec(1));  % Normalize
end

function [x, y] = plot_histogram(vec, freq)
% Outputs x and y so that plot(x, y) generates a normalized histogram
x = zeros(2 * length(freq), 1);
y = zeros(2 * length(freq), 1);
x(1) = vec(1);
x(end) = vec(end);
for i = 2 : length(vec) - 1,
    x(2*i-1) = vec(i);
    x(2*i-2) = vec(i);
end
for i = 1 : length(freq),
    y(2*i) = freq(i);
    y(2*i-1) = freq(i);
end
end
