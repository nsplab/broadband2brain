% Plots data before and after various types of conditioning
% Plots a segment near the beginning and a segment near the end

close all;

% Parameters
data_sets = [1, 2];
elec = {0, 0};  % 0 for all
duration = 0.1;  % Length of each data segment (seconds)
offset = 10;  % Offset from beginning/end (seconds)
buf = 0.01;  % Buffer between data segments (on plot)
cond = {'raw', 'ih'};  % Conditioning types of plot
leg = {'raw', 'highpassed', 'real spikes'};  % Legend
y_sep = 500;  % Separation (y direction) of plots (different conditioning types) on graph
col = {'c', 'b', 'r', 'g', 'k'};
spk_col = {'k', 'b', 'r'};
spk_y = 0;

for data_set = data_sets

    [seg1start, seg2end] = get_total_time(data_set);
    seg1start = seg1start + offset;
    seg2end = seg2end - offset;
    seg1end = seg1start + duration;
    seg2start = seg2end - duration;
    if elec{data_set} == 0
        elec{data_set} = channels_to_use(data_set);
    end

    for i = 1 : length(elec{data_set})
        figure;
        hold on;
        h = [];
        for j = 1 : length(cond)
            [time1 dat1] = get_data(data_set, [elec{data_set}(i)], seg1start, seg1end, cond{j});
            [time2 dat2] = get_data(data_set, [elec{data_set}(i)], seg2start, seg2end, cond{j});
            time2_new = time2 - time2(1) + time1(end) + buf;
            dat1 = dat1 + y_sep * (j-1);
            dat2 = dat2 + y_sep * (j-1);
            fig = plot(time1, dat1, col{j});
            plot(time2_new, dat2, col{j});
            h = [h fig];
        end

        % Plot real spikes
        spikes = real_spikes(data_set, elec{data_set}(i), seg1start, seg1end);
        for j = 1 : length(spikes)
            fig = plot(spikes{j}, ones(size(spikes{j})) * spk_y, [spk_col{j} '.']);
            h = [h fig];
        end
        spikes = real_spikes(data_set, elec{data_set}(i), seg2start, seg2end);
        for j = 1 : length(spikes)
            plot(spikes{j} - time2(1) + time1(end) + buf, ones(size(spikes{j})) * spk_y, [spk_col{j} '.']);
        end

        legend(h, cond);
        title(['data set ' int2str(data_set) ', electrode ' int2str(elec{data_set}(i))]);
        xlabel('time (sec)');
    end

end
