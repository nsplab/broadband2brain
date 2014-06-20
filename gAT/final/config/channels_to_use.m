function [chan] = channels_to_use(data_set)
% Returns list of electrodes to be used
% Use only electrodes with neurons on them

if data_set == 1
    chan = [3, 8, 9, 10, 11, 12, 13, 14, 15, 16];
elseif data_set == 2
    chan = [1, 2, 4, 5, 12, 13, 14];
else
    chan = [5];
end

end
