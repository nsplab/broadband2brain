% Check to see if intervals with spikes have significantly different integrals from ones without

close all;
clear all;

% Parameters
data_set = 1;
electrode = 9;
interval = 0.003;  % Bin width in seconds
threshold = 125;  % For threshold spike detection
delta_th = 100;  % Refractory period measured in indices (for threshold)
cond = 'raw';
%int_handle = @integrate_diode_handle;  % TODO: choose @integrate_handle or @integrate_diode_handle

% Data Segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start;
end_time = data_start + 10;  % Make duration a multiple of interval

corr_vec = [];

for data_set = [1, 2]
    for electrode = channels_to_use(data_set)

        [time, dat] = get_data(data_set, [electrode], start_time, end_time, cond);
        dat_raw = -dat;

        % Highpass  TODO: only use this if processing data directly (without iterative intervals)
        [args, T, window, hp_handle, diode_handle, t0] = uni_params('ker', 1);

        % FILTER SETTINGS
        [A, B] = butter_hp(1000, 2);
        %[A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Causal version
        %hp_handle = @filtfilt_generic;  % Noncausal
        % END

        dat = -hp_handle(dat, args);
        dat_hp = dat;
        
        % Choose 1:
        dat = max(dat, 0);  % Diode (half-wave rectification)
        %dat = abs(dat);  % Full-wave rectification

        % Get spikes times  TODO: choose one of 2 options below
        % Result is thres_t, vector of spike times

        % Option 1: 'real' spikes
        thres_t = real_spikes(data_set, electrode, start_time, end_time, 1);

        % Option 2: run thresholding
        %thres_t = threshold_detect(dat, threshold, delta_th, time);  % TODO: adapt to new

        % Integrate intervals  TODO: choose one of 2 options below
        % Result is integrals(i) = integral over [ endpoints(i) , endpoints(i+1) )

        % Option 1: directly (allows for using highpassed data)

        dt = (time(end)-time(1))/(length(time)-1);
        endpoints = start_time : interval : end_time;
        integrals = zeros(length(endpoints)-1, 1);  
        for i = 1 : length(integrals),
            %int_ind = find(time >= endpoints(i) & time < endpoints(i+1));
            int_ind = find_fast(time, endpoints(i), endpoints(i+1), get_dt());
            integrals(i) = sum(dat(int_ind)) * dt;
        end

        % Option 2: iterate_intervals
        %[t1, result] = iterate_intervals(data_set, [electrode], [start_time], [end_time], cond, start_time, interval, 0, int_handle, []);
        %endpoints = [t1{1} end_time];
        %integrals = result{1};

        % Get spike count in each bin
        spike_vec = bin_spikes(endpoints, thres_t);

        % Vectors to plot for histogram
        x_hist = zeros(1, length(integrals)*2);
        x_hist(1) = endpoints(1);
        x_hist(end) = endpoints(end);
        for i = 2 : length(endpoints)-1,
            x_hist(2*i-2) = endpoints(i);
            x_hist(2*i-1) = endpoints(i);
        end
        y_hist = zeros(1, length(integrals)*2);
        y_count = zeros(1, length(integrals)*2);
        for i = 1 : length(integrals),
            y_hist(2*i-1) = integrals(i);
            y_hist(2*i) = integrals(i);
            y_count(2*i-1) = spike_vec(i);
            y_count(2*i) = spike_vec(i);
        end

        % Plot results
        %{
        figure;
        hold on;
        plot(time, dat, 'c');  % Plot data
        stem(thres_t, 50*ones(size(thres_t)), 'r');  % Plot spikes
        plot(x_hist, y_hist*1000, 'b');  % Plot histogram of integrals
        plot(x_hist, y_count*100, 'g');  % Plot count of spikes
        xlabel('time (sec)');
        legend('data', 'spikes', 'integrals', 'spike count');
        %}

        % Find mean integral for each spike count
        N = max(spike_vec);
        mean_int = zeros(1, N+1);
        n_int = zeros(1, N+1);
        for i = 1 : length(spike_vec),
            mean_int(spike_vec(i)+1) = mean_int(spike_vec(i)+1) + integrals(i);
            n_int(spike_vec(i)+1) = n_int(spike_vec(i)+1) + 1;
        end
        mean_int = mean_int ./ n_int;

        % Plot count vs integral scatterplot
        %{
        figure;
        hold on;
        plot(spike_vec, integrals, 'b.');
        plot(0:N, mean_int, 'r.');
        xlabel('number of spikes');
        ylabel('integral');
        legend('bins', 'averages');
        %}

        % Plot histogram of integral distributions
        fig = figure;
        subplot(1, 2, 1);
        hold on;
        colors = ['b', 'r', 'g'];
        low = min(integrals);
        high = max(integrals)+0.001;
        bins = 50;
        for n = 0 : 1; %N,
            [x, y freq] = make_histogram(low, high, bins, integrals(find(spike_vec == n)));
            plot(x, y, colors(n+1));
            eval(['freq' int2str(n) ' = freq;']);
        end
        labels = {};
        for i = 0 : N,
            labels{i+1} = strcat(int2str(i), ' spikes');
        end

        corr = dot(freq0, freq1);  % Correlation measure
        corr_vec = [corr_vec corr];

        legend(labels);
        xlabel('integral');
        ylabel('frequency');
        title(['data set: ' int2str(data_set) ', elec: ' int2str(electrode) ', corr: ' num2str(corr)]);  % TODO: normalize corr
        axis tight;

        % Plot spike train
        subplot(1, 2, 2);
        hold on;
        y_space = 200;
        plot_time = [thres_t(2) - 0.003, thres_t(2) + 0.003];
        plot_ind = find(time >= plot_time(1) & time < plot_time(2));
        plot(time(plot_ind), dat_raw(plot_ind) - mean(dat_raw(plot_ind)), 'c');
        plot(time(plot_ind), y_space + dat_hp(plot_ind), 'b');
        plot(time(plot_ind), 2*y_space + dat(plot_ind), 'r');
        plot([thres_t(2) thres_t(2)], ylim, 'g--');
        axis tight;

        %saveas(fig, ['~/Desktop/integrals/integrals_' int2str(data_set) '_' int2str(electrode) '.pdf']);

    end
end

av_corr = mean(corr_vec)
