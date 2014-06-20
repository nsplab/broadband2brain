function t_new = impose_delta(t, delta)
% Takes a spike train t (sorted)
% Imposes a refractory period delta
% In the output spike train t_new,
% any 2 spikes must be separated by at least delta

l_new = 1;
t_new = zeros(size(t));
t_new(1) = t(1);
i = 1;
j = 2;
while i < length(t)
    while j <= length(t) && t(j) - t(i) < delta
        j = j + 1;
    end
    if j > length(t)
        break;
    end
    l_new = l_new + 1;
    t_new(l_new) = t(j);
    i = j;
end

t_new = t_new(1 : l_new);

end