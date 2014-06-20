function [t1, result] = iterate_intervals(data_set, electrodes, start_times, end_times, cond, t0, T, window, handle, args)
% Generic function for dealing with large portions of data piece by piece
% Takes in 1 or more arbitrary data segments (in terms of time)
% Divides the segments into intervals of length T
% For each interval, passes data to the handler function handle
% Concatenates results for each interval
% See the folder /handle for possible handler functions

% Inputs:
% data_set - data set number (1 or 2)
% electrodes - vector of electrode numbers to process
% start_times - vector of segment start times (sorted)
% end_times - vector of segment end times
% cond - whether data should be conditioned (see config/data_file.m)
% t0 - offset, i.e. first interval of size T starts at t0
% T - interval size
% window - extra data padding required on either size of interval (in seconds)
% handle - handler function passed as @function_name
% args - struct of extra parameters for handler function

% Output to handler function (for each interval):
% data_set
% electrode
% seg_start - segment start time
% seg_end - segment end time
% t1 - interval start time
% T - length of interval
% time - time vector
% dat - data vector
% args - struct of extra (handler-specific) parameters

% Return values of handler function (for each interval):
% res - result, matrix (should be a fixed height for a given handler)

% Return values:
% t1{segment}(:) = list of interval start times in this segment
% result - concatenation of results from intervals
% * result{electrode_num, segment_num} = horizontal concatenation of (matrix) results from intervals in segment segment_num

assert(length(start_times) == length(end_times), 'Error: seg_start and seg_end must have the same length');

num_segs = length(start_times);
num_elecs = length(electrodes);
cur_file_num = 0;  % Most recently opened file num
time_vec = [];  % time_vec(elec) = current time vector
data_vec = [];  % data_vec(elec_num, :)
t1 = cell(num_segs, 1);
result = cell(num_elecs, num_segs);

data_to_pass = cell(num_elecs, 1);

% for concatenating results into large matrix
res_length = zeros(length(electrodes));  % indices allocated
res_used = zeros(length(electrodes));  % indices used

for cur_seg = 1 : num_segs
    seg_start = start_times(cur_seg);
    seg_end = end_times(cur_seg);
    seg_start_file = seg_lookup(data_set, seg_start - window);  % File num where this segment starts

    % Make sure correct file is loaded for start of segment
    if cur_file_num ~= seg_start_file
        % Load seg_start_file
        [time_vec, data_vec] = load_data_file(data_set, electrodes, seg_start_file, cond);
        cur_file_num = seg_start_file;
    end

    % Partition segment into intervals
    first_int_start = seg_start - mod(seg_start - t0, T);
    int_start = first_int_start : T : seg_end;  % Contains interval start times
    if int_start(end) == seg_end
        int_start = int_start(1 : end-1);
    end
    t1{cur_seg} = int_start;

    print_counter = 0;

    % Process each interval
    for i = 1 : length(int_start)
        
        ticIDiter = tic;

        print_counter = print_counter + 1;
        if print_counter == 100
            disp([int2str(i) ' of ' int2str(length(int_start)) ', ' num2str(i/length(int_start)*100) '%   iter time = ' num2str(iter_time) '   concat time = ' num2str(concat_time) '   corr time = ' num2str(iter_time - concat_time)]);
            print_counter = 0;
        end

        %disp('start');

        %ticID = tic;

        % If interval goes beyond currently loaded data, load more data
        while int_start(i) + T + window > time_vec(end) && cur_file_num+1 <= num_data_segs(data_set)
            % Decide how much of previous data needs to be kept
            keep_ind = find(time_vec >= int_start(i) - window);
            % Load next file
            [time_new, data_new] = load_data_file(data_set, electrodes, cur_file_num + 1, cond);
            time_vec = [time_vec(keep_ind) time_new];
            data_vec = [data_vec(:, keep_ind) data_new];
            cur_file_num = cur_file_num + 1;
        end

        %toc(ticID)
        %ticID = tic;

        % Isolate data for this interval
        % Result: int_ind is set to vector of indices in time_vec in range [int_start(i) - window, int_start(i) + T + window)
        %int_ind = find(time_vec >= int_start(i) - window & time_vec < int_start(i) + T + window);  % Original - slow
        int_ind = find_fast(time_vec, int_start(i) - window, int_start(i) + T + window, get_dt());  % Optimized

        %toc(ticID)
        %ticID = tic;

        % Make sure this interval actually contains some relevant data
        % Need this line in case interval starts right at end of segment
        if time_vec(int_ind(1)) >= seg_end
            continue;
        end

        % Call handler
        for elec_num = 1 : length(electrodes)

            ticID2 = tic;

            args.data_set = data_set;
            args.channel = electrodes(elec_num);  % pass channel num
            [res next_data] = handle(data_set, electrodes(elec_num), seg_start, seg_end, int_start(i), T, time_vec(int_ind), data_vec(elec_num, int_ind), args, data_to_pass{elec_num});

            %run_handler_time = toc(ticID2)
            
            data_to_pass{elec_num} = next_data;  % Pass on data

            data_pass_time = toc(ticID2);

            if i == 1
                % preallocate memory, guess how big array is likely to
                % become
                result{elec_num, cur_seg} = zeros(size(res, 1), length(int_start)*2);
                res_length(elec_num) = length(int_start)*2;
                res_used(elec_num) = 0;
                %result{elec_num, cur_seg} = [result{elec_num, cur_seg} res];
            end
            
            if res_used(elec_num) + size(res, 2) >= res_length(elec_num)
                % need to resize array, double its length
                result{elec_num, cur_seg} = [result{elec_num, cur_seg} zeros(size(result{elec_num, cur_seg}))];
                res_length(elec_num) = res_length(elec_num) * 2;
            end
            
            base_ind = res_used(elec_num) + 1;
            result{elec_num, cur_seg}(:, base_ind : base_ind + size(res, 2) - 1) = res;
            res_used(elec_num) = res_used(elec_num) + size(res, 2);
            
            if i == length(int_start)
                % done, trim result
                result{elec_num, cur_seg} = result{elec_num, cur_seg}(:, 1 : res_used(elec_num));
            end
                
            total_time = toc(ticID2);
            concat_time = total_time - data_pass_time;

        end

        iter_time = toc(ticIDiter);

    end

end

end
