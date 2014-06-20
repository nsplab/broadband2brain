function [y, sample_times] = sample_conv(time, x, t0, T, L, h, args)
% Takes L samples in interval [t0, T] by convolution with h(t)
% Half as much spacing at ends as between 2 samples

T_s = T / L;
sample_times = linspace(t0 + T_s/2, t0 + T - T_s/2, L);

y = convolve(time, sample_times, x, h, args);

end
