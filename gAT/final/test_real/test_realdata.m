function [error runtime] = test_realdata(method_name, paramID)
% Runs a method on real data
% Computes reconstruction error (samples vs model-reconstructed samples)
% Reports runtime
% If debug = 1, plots data, reconstructed spikes, samples

debug = 1;

% Parameters
data_set = 1;
electrode = 9; %channels_to_use(data_set);
cond = 'raw';
duration = 1;

% Find data segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start;
end_time = start_time + duration;

% Get params from paramID
[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);

% Clear file for saved samples
filename = [root_dir() 'mat_files/samples.mat'];
t_samples = [];
samples = [];
samples_model = [];
runtime = [];
save(filename, 't_samples', 'samples', 'samples_model', 'runtime');

% Reconstruct spikes
[t1, result] = run_generic(data_set, [electrode], [start_time], [end_time], hp_handle, diode_handle, method_name, args, start_time, T, window);
[t, c, t1, sigma] = extract_result_single(t1, result);
t = t{1};
c = c{1};
t1 = [t1 t1(end)+T];
sigma = sigma{1};

% Real spikes
spikes = real_spikes(data_set, electrode, start_time, end_time, 1);

% Samples
load(filename);
error = rec_error(samples, samples_model);  % Reconstruction error
runtime = mean(runtime);

% Plot results
if debug == 1
    sc = 1000;
    fig = figure;
    hold on;
    [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
    dat = -hp_handle(dat, args);  % Highpass
    plot(time, dat, 'c');  % Real data
    stem(t, c*sc, 'k');  % Reconstructed spikes
    plot(t1, ones(size(t1))*-25, 'b.');  % Sampling times
    plot(spikes, ones(size(spikes))*-15, 'k.');  % Real spikes
    stem(t, sigma, 'm');  % Plot sigma
    plot(t_samples, samples, 'm.-');  % Samples
    plot(t_samples, samples_model, 'b.-');  % Model samples

    xlim([time(1), time(end)]);

    leg = legend('data', 'reconstructed spikes', 'sampling times', 'real spikes', 'sigma', 'samples', 'model samples');
    set(leg, 'Location', 'SouthOutside');
    xlabel('time (sec)');
end

end