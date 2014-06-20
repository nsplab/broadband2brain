function [setup_obj] = setup_gibbs(T, window, args)
% Precomputes values for gibb's sampling method

L = args.L;

% Set up for convolution
convolve_setup_2(T, window, L, @gaussian, args);

setup_obj = struct();

end
