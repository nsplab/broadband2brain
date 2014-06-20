function h = gaussian(t, args)
% args needs args.sigmah - kernel width

h = exp(-1/(2*args.sigmah^2)*t^2);

end
