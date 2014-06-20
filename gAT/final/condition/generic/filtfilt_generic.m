function x_out = filtfilt_generic(dat, args)

x_out = filtfilt(args.B_filter, args.A_filter, dat);

end
