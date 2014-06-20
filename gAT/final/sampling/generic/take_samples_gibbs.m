function [y sample_times] = take_samples_gibbs(time, dat, t1, T, args)
% Takes samples by convolution with gaussian kernel
[y, sample_times] = convolve_fast(dat, t1);

end
