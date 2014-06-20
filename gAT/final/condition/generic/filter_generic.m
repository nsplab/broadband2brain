function x_out = filter_generic(dat, args)

x_out = filter(args.B_filter, args.A_filter, dat);

end
