function ind = find_fast(time_vec, low_target, high_target, dt)
% Returns indices of (sorted) vector vec in range [low_target, high_target)
% Alternative to matlab's find function

low_guess = 1 + round((low_target - time_vec(1)) / dt);
if low_guess < 1
    low_guess = 1;
end
if low_guess > length(time_vec)
    low_guess = length(time_vec);
end
while low_guess >= 1 && time_vec(low_guess) >= low_target
    low_guess = low_guess - 1;
end
low_guess = low_guess + 1;
while low_guess <= length(time_vec)-1 && time_vec(low_guess) < low_target
    low_guess = low_guess + 1;
end

high_guess = 1 + round((high_target - time_vec(1)) / dt);
if high_guess < 1
    high_guess = 1;
end
if high_guess > length(time_vec)
    high_guess = length(time_vec);
end
while high_guess >= 2 && time_vec(high_guess) >= high_target
    high_guess = high_guess - 1;
end
while high_guess <= length(time_vec) && time_vec(high_guess) < high_target
    high_guess = high_guess + 1;
end
high_guess = high_guess - 1;

ind = low_guess : high_guess;

end
