function [y xs] = succ_int_2(x, dt, L)
% Takes L successive integral samples
% Outputs full signal for each integral (not just samples at endpoint)

N = size(x, 2);

xs = cell(1, L);
x_si = zeros(L, N);

for i = 1:N
    x_si(:, i) = succ_int(x(1:i), dt, L);
end
y = x_si(:, N);

for i = 1:L
    xs{i} = x_si(i, :);
end

%xs = cell(1, L);
%
%y = zeros(L, 1);
%x_int = integrate(x, dt);
%xs{1} = x_int;
%y(1) = x_int(end);
%for i = 2 : length(y),
%    x_int = integrate(x_int, dt);
%    xs{i} = x_int;
%    y(i) = x_int(end);
%end

end
