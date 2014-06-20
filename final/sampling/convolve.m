function [y] = convolve(time, sample_times, x, h, args)
% Convolves x with h(t) and samples at points in time given in sample_times
% h must accept a single array of the form [t, arg1, arg2, ...]
% Pass h as function handle, e.g. @normalized_sinc

% NOTE: this is the old convolve function that runs very slowly

%ticID_convolve = tic;

y = zeros(length(sample_times), 1);
dt = (time(end) - time(1)) / (length(time) - 1);
for j = 1 : length(y)
    d = 0;
    for k = 1 : length(x)
        d = d + x(k) * h(sample_times(j) - time(k), args) * dt;
    end
    y(j) = d;
end

%time_convolve = toc(ticID_convolve)

end
