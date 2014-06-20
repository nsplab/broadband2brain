% script for running LFP analysis on all channels
% run_LFP_2 saves results in .mat files

close all;

data_set = 1;
compute = 1;

channels = channels_to_use(data_set);

for i = 1 : length(channels)
    run_LFP_2(data_set, channels(i), compute);
end