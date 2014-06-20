function y = succ_int(x, dt, L)
% Takes L successive integral samples

%N = size(x, 2);
%T = N * dt;
%t = dt * ((1:N) - 0.5), L;
%y = zeros(L, 1);
%for i = 1:L
%    for j = 1:N
%        y(i) = y(i) + 

M = size(x, 1);
N = size(x, 2);
A = repmat((1:L)', 1, N);
y = zeros(L, M);
T = N * dt;
t1 = repmat(dt * (0:(N-1)), L, 1);
t2 = repmat(dt * (1:N), L, 1);
C = ((T - t1) .^ A - (T - t2) .^ A) ./ factorial(A);
for i = 1:M
    B = repmat(x(i, :), L, 1);
    y(:, i) = sum(B .* C, 2);
end


%N = size(x, 2);
%A = repmat((0:(L-1))', size(x));
%B = repmat(x, L, 1);
%T = N * dt;
%t = repmat(dt * ((1:N) - 0.5), L, 1);
%y = dt * sum(B .* (T - t) .^ A ./ factorial(A), 2);


%y = zeros(L, 1);
%x_int = integrate(x, dt);
%y(1) = x_int(end);
%for i = 2 : length(y),

%    x_int = integrate(x_int, dt);
%    y(i) = x_int(end);
%end

end
