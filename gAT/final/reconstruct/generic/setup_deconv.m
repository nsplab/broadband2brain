function [setup_obj] = setup_deconv(T, window, args)
% Precomputes values for generic kernel
% args.h is function handler for kernel
% args.h_args is args for kernel

% copy/pasted from setup_ker.m
% not all parts of this are actually used - this is kind of a hack

L = args.L;
h = args.h;

% Set up for convolution
convolve_setup(T, window, L, h, args);

T_s = T / L;  % Spacing of samples

setup_obj = struct();

end