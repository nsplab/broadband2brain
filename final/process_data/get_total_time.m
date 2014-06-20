function [start_time, end_time, ff_start, ff_end, wash_start, wash_end] = get_total_time(data_set)
% Returns start and end time of data
% If relevant (nemo channel 5 analysis), returns start/end for
% baseline/ff/washout
% Run config_lookup.m before this will work

filename = [root_dir() 'process_data/data_lookup.mat'];  % Must be the same as in config_lookup.m
load(filename);

start_time = start_time(data_set);
%end_time = end_time(data_set);

% In case unused
ff_start = 0;
ff_end = 0;
wash_start = 0;
wash_end = 0;

% End time at start of force field data
if data_set == 1
    end_time = end_time(data_set);  % All data (control day)
%elseif data_set == 2
%    end_time = 1039.1;  % Baseline (force field day)
%{
elseif data_set == 3
    end_time = 1245.0;
elseif data_set == 4
    end_time = 1056.0;
elseif data_set == 5
    end_time = 1262.5;
else
%}
else
    % nemo channel 5 analysis
    [end_time, ff_start, ff_end, wash_start, wash_end] = test_baseline_trials(data_set);
end

end
