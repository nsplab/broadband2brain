%  spike_example(1, [8 9 13 16], 0.2095, 0.07, [0.0281 0.0084, 0.0378 0.0683]-0.001, 0.002)
function [] = spike_example(data_set, electrodes, start_offset, T, zoom_start, T_zoom)

% Data segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start + start_offset;
end_time = start_time + T;

N = length(electrodes);

% Get data
for i = 1:N

    [time, dat] = get_data(data_set, electrodes(i), start_time, end_time, 'raw');
    time = time - start_time;  % start first interval at time 0
    
    % Highpass
    [A, B] = butter_bp(300, 6000, 2);  % Change filter order?
    args.A_filter = A;
    args.B_filter = B;
    hp_handle = @filter_generic;  % Use causal version
    dat = hp_handle(dat, args);
    
    % Raw data (AT)
    subplot(N + 1, N, (1:N) + N * (i - 1));
    hold on;
    plot(time*1000, -dat, 'Color', [21 51 173]/255);

    spikes = real_spikes(data_set, electrodes(i), start_time, end_time, 0);
    %yl = [min(-dat) max(-dat)];
    frac = 0.05;
    range = max(-dat) - min(-dat);
    yl = [min(-dat)-frac*range max(-dat)+frac*range];
    xl = xlim;
    spikes = spikes{1};
    spike_offset = spikes - start_time;
    i
    spike_offset
    for j = 1:length(spikes)
        t = spikes(j) - start_time;
        scatter([t] * 1000, yl(1), 'black');
    end
    xlim(xl);
    ylim(yl);
    box on;
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);



    if nargin == 6
	valid = ((time >= zoom_start(i)) & (time <= (zoom_start(i) + T_zoom)));
	plot([min(time(valid)) min(time(valid))]*1000, yl, 'black');
	plot([max(time(valid)) max(time(valid))]*1000, yl, 'black');

        subplot(N + 1, N, N * N + i);
	hold on;
	plot(time(valid)*1000, -dat(valid), 'Color', [21 51 173]/255);
	scatter(1000*(zoom_start(i)+T_zoom/2), yl(1), 'black');
	xlim([min(time(valid)) max(time(valid))]*1000);
        ylim(yl);
	box on;
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end

end

