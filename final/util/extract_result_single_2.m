function [t, c, t1, sigma] = extract_result_single_2(t1, result, elecs)
% Extracts result produced by iterate_intervals for a single segment of data
% Retains channel number info

t1 = t1{1};

% Extract t, c from result
t = cell(size(max(elecs), 1), 1);
c = cell(size(max(elecs), 1), 1);
sigma = cell(size(max(elecs), 1), 1);
for i = 1 : length(result)
    t{elecs(i)} = result{i}(1, :);
    c{elecs(i)} = result{i}(2, :);
    sigma{elecs(i)} = result{i}(3, :);
end

end
