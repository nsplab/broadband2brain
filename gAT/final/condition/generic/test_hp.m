function [A, B] = test_hp(high_cutoff, N)
% Returns parameters A, B

% Inputs:
% high_cutoff - cutoff freq in Hz
% N - order of filter

% Parameters
sampl_freq = 1 / get_dt();  % Sampling frequency in Hz

Wn = high_cutoff / (sampl_freq / 2);  % Cutoff freq normalized to half sample rate

[B, A] = cheby1(N, Wn, 'high');

end
