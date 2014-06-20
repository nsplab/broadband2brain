function [t, c, t1, sigma] = extract_result_single(t1, result)
% Extracts result produced by iterate_intervals for a single segment of data

t1 = t1{1};

% Extract t, c from result
t = cell(size(result, 1), 1);
c = cell(size(result, 1), 1);
sigma = cell(size(result, 1), 1);
for i = 1 : length(result)
    t{i} = result{i}(1, :);
    c{i} = result{i}(2, :);
    sigma{i} = result{i}(3, :);
end

end
