function [A, B] = butter_hp(high_cutoff, N)
% Returns parameters A, B for butterworth highpass filter

% Inputs:
% high_cutoff - cutoff freq in Hz
% N - order of filter

% Parameters
sampl_freq = 1 / get_dt();  % Sampling frequency in Hz

Wn = high_cutoff / (sampl_freq / 2);  % Cutoff freq normalized to half sample rate

[B, A] = butter(N, Wn, 'high');

end
