function [y sample_times] = take_samples_RMSE(time, dat, t1, T, args)
% Takes samples by convolution with generic kernel
% same as take_samples_ker.m

[y, sample_times] = convolve_fast(dat, t1);

end
