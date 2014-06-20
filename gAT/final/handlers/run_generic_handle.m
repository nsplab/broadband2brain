function [res next_data] = run_generic_handle(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data)
% Runs generic reconstruction method and returns spike times/amplitudes
% res(1, :) = times
% res(2, :) = amplitudes
% res(3, 1) = sigma
% next_data is a struct to be passed on to next interval (containing spike t,c in this interval)

ticID1 = tic;

debug = 0;  % Set to 1 to make plots

% Apply flip
dat = -dat;

% Apply highpass and diode
hp_handle = args.hp_handle;
diode_handle = args.diode_handle;
if ~isnumeric(hp_handle)
    dat = hp_handle(dat, args);  % Highpass
end
if ~isnumeric(diode_handle)
    dat = diode_handle(dat, args);  % Diode
end

% TESTING: plot
if debug
    test_c_stability(data_set, elec, seg_start, seg_end, t1, T, time, dat, args, prev_data);
end

%ticID_s = tic;

% Take samples
sampleID = tic;
sample_func = eval(['@take_samples_' args.method_name]);
[y, sample_times] = sample_func(time, dat, t1, T, args);
%sample_time = toc(sampleID)
%[y] = sample_func(time, dat, t1, T, args);

%samp_time = toc(ticID_s)
%ticID_r = tic;

% Run reconstruction
rec_func = eval(['@reconstruct_' args.method_name]);
[t, c, sigma, elapsed_time] = rec_func(y, T, args, prev_data);

%elapsed_time

%rec_time = toc(ticID_r)

% TESTING: run test_likelihood
if debug
    [t_vec, c_vec, sig_vec, like_vec] = test_likelihood(y, T, args, prev_data);
    t_vec = t_vec + t1;
end

% Pass spikes on to next interval
next_data = struct();
next_data.t = t;
next_data.c = c;

t = t + t1;

if length(t) == 0
    t = zeros(0, 1);
end

% Limit spike times to [seg_start, seg_end]
ind = find(t >= seg_start & t < seg_end);
t = t(ind, 1);
c = c(ind, 1);
    
% Output
res = zeros(3, length(t));

res(1, :) = t';
res(2, :) = c';
res(3, :) = sigma * ones(size(res(3, :)));

% TESTING: plot

%{

if length(c) == args.K  % at very end of interval c will have length 1, ignore

    % get model-reconstructed samples
    y_model = model_samples(t1, sample_times, T, args, prev_data, t, c, sigma);
    
    % Save samples to file for plotting later
    conID = tic;
    filename = [root_dir() 'mat_files/samples.mat'];
    load(filename);
    t_samples = [t_samples sample_times];
    samples = [samples ; y];
    samples_model = [samples_model ; y_model];
    runtime = [runtime elapsed_time];
    handle_concat_time = toc(conID)
    saveID = tic;
    save(filename, 't_samples', 'samples', 'samples_model', 'runtime');
    save_time = toc(saveID)
    
    if debug
        sinc_t = time(1) : 0.00001 : time(end);
        sinc_x = zeros(size(sinc_t));
        for i = 1 : length(sinc_x)
            sinc_x(i) = args.h(sinc_t(i) - mean(sinc_t), args);
        end
        figure;
        subplot(2, 1, 1); % reconstructed spikes, samples, kernel
        hold on;
        plot(time, dat, 'c');
        if ~isnan(sample_times)
            plot(sample_times, y*5, 'm.-'); % plot samples
        end

        % plot model-reconstructed samples
        plot(sample_times, y_model*5, 'b.-');

        stem(t, c*10000, 'k');
        plot([t1, t1+T], [-5, -5], 'b.');
        y_lim = ylim;
        ker_scale = y_lim(2) / max(sinc_x);
        plot(sinc_t, sinc_x*ker_scale/2, 'r');
        legend('data', 'samples', 'model samples', 'reconstructed', 'interval', 'kernel');
        xlabel('time (sec)');
        title(['sigma = ' num2str(sigma)]);
        xlim([t1, t1+T]);
        ylim([-20, 240]);

        subplot(2, 1, 2); % c, t, L as a function of spike position
        hold on;
        plot(t_vec, c_vec, 'k.-');
        plot(t_vec, sig_vec * max(c_vec) / max(sig_vec), 'm.-');
        plot(t_vec, like_vec * max(c_vec) / max(like_vec), 'b.-');
        legend('c', 'sigma', 'likelihood');
        axis tight;
        ylim([-0.05, 0.05]);
    end
end

%}

%handle_time = toc(ticID1)

end
