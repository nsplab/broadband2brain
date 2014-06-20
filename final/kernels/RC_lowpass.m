function [h] = RC_lowpass(t, args)
% Impulse response of a first-order RC lowpass filter
% From pg 240 eq (3.146) of Signals and Systems by Oppenheim and Willsky
% args = [RC]

%RC = args(1);
RC = args.RC;

if t < 0,
    h = 0;
else
    h = 1/RC * exp(-t / RC);
end

end
