function [] = reconstruct_and_save(method_name, paramID, data_sets, alldata_flag)
% Runs reconstruction method on entire data set and saves spike times/amplitudes to a .mat file
% alldata_flag is optional - set to 1 to run on all (ff, wash) data

if nargin < 4
    alldata_flag = 0
end

disp('NEW RUN');
method_name
paramID

[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);

for data_set = data_sets

    data_set

    [start_time, end_time] = get_total_time(data_set);
    if alldata_flag == 1
        [start_time, base_end, ff_start, ff_end, wash_start, end_time] = get_total_time(data_set);
    end

    start_time = start_time + 1;  % Bug fix

    electrodes = channels_to_use(data_set);
    %electrodes = [16];  % TESTING
    
    [t1, result] = run_generic(data_set, electrodes, start_time, end_time, hp_handle, diode_handle, method_name, args, t0, T, window);
    
    [t, c, t1, sigma] = extract_result_single_2(t1, result, electrodes);

    filename = [root_dir() 'reconstruct/saved_spikes/spk_' method_name '_' int2str(data_set) '_' int2str(paramID) '.mat'];
    save(filename, 't', 'c', 'sigma');  % TODO: add t1?

end

end