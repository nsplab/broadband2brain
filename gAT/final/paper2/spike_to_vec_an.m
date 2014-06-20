function vec = spike_to_vec_an(spike_times, spike_amplitudes, T, bin_size, start_time, num_bins)
% turns spike times into a spike vector
% 1 reconstructed spike per T-period
% bin_size (~ 1 ms) is for answer
% T must be a multiple of bin_size
% start_time is the start of a bin, not necessarily the start of a T-period

end_time = start_time + bin_size*(num_bins-1/2);

% make start_time the start of a T-period
leading_bins = round(mod(start_time, T) / bin_size);
start_time = start_time - leading_bins*bin_size;

start_times_T = start_time : T : end_time;
multiple = round(T / bin_size);

vec = zeros(1, num_bins);

% find index of spike_times corresponding to first T-period
av_start_time = mean(spike_times' - (0:T:((length(spike_times)-1/2)*T)));
start_ind = round((start_time - av_start_time)/T + 1/2);  % should be +1 but corrected below

for i = 1 : length(start_times_T)
    if spike_amplitudes(start_ind + i) > 0
        t = spike_times(start_ind + i) - start_times_T(i);
        t_ind = round(t / bin_size + 1/2);
        if t_ind < 1
            t_ind = 1;
        elseif t_ind > multiple
            t_ind = multiple;
        end
        ind = t_ind + (i-1)*multiple - leading_bins;
        if ind >= 1 && ind <= num_bins
            vec(ind) = 1;
        end
    end
end

% enforce 1 ms refractory period (1 bin)
for i = 2 : num_bins
    if vec(i) == 1 && vec(i-1) == 1
        vec(i) = 0;
    end
end

end