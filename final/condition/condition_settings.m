function [] = condition_settings(data_sets, cond)
% Settings for different types of data conditioning are here
% Call this function in order to condition and re-save data

% condition_data(data_sets, highpass_func, highpass_args, diode_func, diode_args, cond, window, flip)

if strcmp(cond, 'raw')  % Raw data
    % Do nothing because raw data is already there
elseif strcmp(cond, 'ih')  % Ideal highpass (noncausal butterworth)
    condition_data(data_sets, @butter_noncausal, [], 0, [], cond);
elseif strcmp(cond, 'ihd')  % Ideal highpass and diode
    condition_data(data_sets, @butter_noncausal, [], @ideal_diode, [], cond);
elseif strcmp(cond, 'ch')  % Causal highpass
    condition_data(data_sets, @butter_causal, [], 0, [], cond);
elseif strcmp(cond, 'ch1')  % 1st order causal butterworth
    condition_data(data_sets, @butter_causal, [-1, -1, 1], 0, [], cond);
end

end
