function [x_hp] = butter_noncausal(x, args)
% High-pass the data using filtfilt (noncausal)
% args = [sampl_freq, high_cutoff, order]

% Parameters
if length(args) < 1 || args(1) == -1
    sampl_freq = 1/get_dt();  % Sampling frequency in Hz
else
    sampl_freq = args(1);
end
if length(args) < 2 || args(2) == -1
    high_cutoff = 400;  % High-pass cutoff in Hz
else
    high_cutoff = args(2);
end
if length(args) < 3 || args(3) == -1
    N = 10;  % Order of filter
else
    N = args(3);
end

Wn = high_cutoff / (sampl_freq / 2);  % Cutoff freq normalized to half sample rate

[B, A] = butter(N, Wn, 'high');  % Uses Butterworth filter

x_hp = filtfilt(B, A, x);

B
A

end
