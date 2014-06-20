function [t1, result] = run_generic(data_set, electrodes, start_times, end_times, hp_handle, diode_handle, method_name, args, t0, T, window)
% Generic function to run any reconstruction method

% Naming conventions:
% method_name is a string such as 'int' or 'sinc'
% The setup function is called setup_<method_name>
% The reconstruction function is called reconstruct_<method_name>
% The sampling function is called take_samples_<method_name>

% Inputs:
% start_times, end_times - start and end times for segments for iterate intervals
% hp_handle - highpass function handler (0 for none)
% diode_handle - diode function handler (0 for none)
% method_name - string describing method (see naming conventions above)
% args - arguments for reconstruction method, e.g. K, L, delta_t, numiter, t_resolution, kernel, parameters for kernel/highpass/etc
% t0 - start time of first interval (of length T)

% Outputs:
% t - spike times, sorted
% c - amplitudes
% t1 - interval start times
% sigma - sigma values corresponding to t, c above (e.g. contains duplicate entries if K = 2)

% Setup function
setup_func = eval(['@setup_' method_name]);
setup_obj = setup_func(T, window, args);

% Append to args
args.setup_obj = setup_obj;
args.hp_handle = hp_handle;
args.diode_handle = diode_handle;
args.method_name = method_name;

% Run iterate_intervals
[t1, result] = iterate_intervals(data_set, electrodes, start_times, end_times, 'raw', t0, T, window, @run_generic_handle, args);

end
