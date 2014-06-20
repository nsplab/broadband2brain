function [t_new t_ind c_new sigma_new] = select_spikes(t, start_time, end_time, c, sigma)
% Selects spikes in time interval quickly using binary search
% If optional params c, sigma provided, prunes them accordingly


