function [t, c] = sort_t_c(t_in, c_in)
% Re-orders t,c pairs by increasing t

t = t_in;
c = c_in;

for i = 1 : length(t)
    min_t = t(i);
    min_ind = i;
    for j = i + 1 : length(t)
        if t(j) < min_t
            min_t = t(j);
            min_ind = j;
        end
    end
    j = min_ind;
    % Swap indices i and j
    temp_t = t(i);
    temp_c = c(i);
    t(i) = t(j);
    c(i) = c(j);
    t(j) = temp_t;
    c(j) = temp_c;
end

end
