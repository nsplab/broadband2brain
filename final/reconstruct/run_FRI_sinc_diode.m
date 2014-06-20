function [t, c, t1, sigma] = run_FRI_sinc_diode(data_set, electrodes, start_time, end_time, cond, t0, T, L, K, B, delta_t, numiter, t_resolution)
% Runs FRI construction for sinc samples

% Inputs:
% data_set
% electrodes - vector
% start_time, end_time

% Outputs:
% t{elec_num} = [list of times]
% c{elec_num} = [list of amplitudes]
% sigma{elec_num} = [list of sigmas (may include 0's if K > 1)]
% t1 - vector of interval start times

% Parameters
window = 0.01;

% Set up for convolution
convolve_setup(T, window, L, @normalized_sinc, [B]);

[cov_inv, mu, t_sampled, B_full, B_const] = FRI_sinc_diode_setup(L, T, B, t_resolution);

% args is a struct with fields L, K, delta_t, numiter, cov_inv, mu, t_sampled, B_full, B_const
args = struct();
args.L = L;
args.K = K;
args.B = B;
args.delta_t = delta_t;
args.numiter = numiter;
args.cov_inv = cov_inv;
args.mu = mu;
args.t_sampled = t_sampled;
args.B_full = B_full;
args.B_const = B_const;

[t1, result] = iterate_intervals(data_set, electrodes, [start_time], [end_time], cond, t0, T, window, @run_sinc_diode_handle, args);
t1 = t1{1};

% Extract t, c from result
t = cell(length(electrodes), 1);
c = cell(length(electrodes), 1);
sigma = cell(length(electrodes), 1);
for i = 1 : length(result)
    t{i} = result{i}(1, :);
    c{i} = result{i}(2, :);
    sigma{i} = result{i}(3, :);
end

end
