function spikes = real_spikes(data_set, elec, time_start, time_end, combine)
% Returns spike times of 'real' spikes

% Inputs:
% data_set
% elec - single electrode
% time_start and time_end are optional
% combine - 0 or 1 for whether separate neurons should be combined (default 0)

% Outputs:
% if combine == 0, spikes{neuron_num} = [list of spike times]
% if combine == 1, [list of spike times]

if nargin < 5
    combine = 0;
end

[d, info] = data_file(data_set, 0, 0);
load(info);  % Loads spk
spikes = {};

i = 1;
for j = 1 : length(spk)
    if spk(j).elec == elec
        spk_t = spk(j).data;
        if nargin >= 3
            spk_t = spk_t(find(spk_t >= time_start));
        end
        if nargin >= 4
            spk_t = spk_t(find(spk_t < time_end));
        end
        spikes{i} = spk_t;
        i = i + 1;
    end
end

if combine == 1 && length(spikes) == 1
    spikes = spikes{1};
else if combine == 1

    % Start of implementation for more optimized method        
    %{
    d = 0;
    for i = 1 : length(spikes)
        d = d + length(spikes{i});
    end
    spk_comb = zeros(d, 1);
    inds = cell(size(spikes));
    for i = 1 : length(inds);
        inds{i} = 1;
    end
    for res_ind = 1 : d
        
    spikes = spk_comb;
    %}

    spk_comb = [];
    for i = 1 : length(spikes)
        spk_comb = [spk_comb ; spikes{i}];
    end
    spk_comb = sort(spk_comb);
    spikes = spk_comb;

end

end
