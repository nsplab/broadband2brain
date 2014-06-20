function [h] = RLC(t, args)
% Impulse response of series RLC circuit
% Parameters: args.s1 and args.s2 (both negative)
% I.R. = A * exp(s1*t) - A * exp(s2*t) for t >= 0
% A = s1*s2/(s2-s1)

if t < 0,
    h = 0;
else
    A = args.s1*args.s2/(args.s1-args.s2);
    h = A * (exp(args.s1*t) - exp(args.s2*t));
end

end
