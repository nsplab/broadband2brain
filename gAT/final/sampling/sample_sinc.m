function [y sample_times] = sample_sinc(time, x, t0, T, B, L)
% Takes samples by convolution with sinc kernel

[y, sample_times] = sample_conv(time, x, t0, T, L, @normalized_sinc, [B]);

end
