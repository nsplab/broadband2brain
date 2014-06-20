function [y sample_times] = take_samples_ML(time, dat, t1, T, args)
% Takes samples by convolution with generic kernel

[y, sample_times] = convolve_fast(dat, t1);

end
