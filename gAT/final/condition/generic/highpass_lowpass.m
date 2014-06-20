function x_out = highpass_lowpass(dat, args)
% apply highpass and then lowpass
% for analog method, hopefully allows spike sorting based on width

x = filter(args.B_filter, args.A_filter, dat);
x_out = filter(args.B_lpfilter, args.A_lpfilter, x);

end
