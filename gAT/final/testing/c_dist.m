% Plots distribution of reconstructed amplitudes c

close all;
clear all;

method_name = 'int';
paramID = 1;
data_set = 1;
elec = 16;

filename = [root_dir() 'reconstruct/saved_spikes/spk_' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
load(filename);

amp = c{elec};

low = 0;
high = 0.2
bins = 50;
[x, y] = make_histogram(low, high, bins, amp);

figure;
plot(x, y);
xlabel('c');
ylabel('freq');
title(['method: ' method_name ', data set: ' int2str(data_set) ', elec: ' int2str(elec)]);
