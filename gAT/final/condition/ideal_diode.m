function x_new = ideal_diode(x, args)
% Applies an ideal diode to the data, eliminating negative entries
% args = []

x_new = max(x, 0);

end