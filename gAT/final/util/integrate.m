function x_int = integrate(x, dt)
% Integrates a vector, i.e. returns a vector of partial sums (times dt)
x_int = zeros(length(x), 1);
sum = 0;
for i = 1 : length(x),
    sum = sum + x(i) * dt;
    x_int(i) = sum;
end
end
