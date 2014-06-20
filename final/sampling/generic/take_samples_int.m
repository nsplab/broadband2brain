function [y sample_times] = take_samples_int(time, x, t1, T, args)
% Takes L successive integral samples
% Assumes window = 0

dt = get_dt();
L = args.L;

ind = find_fast(time, t1, t1+T, dt);
x = x(ind);

y = zeros(L, 1);
x_int = integrate(x, dt);
y(1) = x_int(end);
for i = 2 : length(y),
    x_int = integrate(x_int, dt);
    y(i) = x_int(end);
end

sample_times = NaN;  % For consistency

end
