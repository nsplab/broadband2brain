function h = normalized_sinc(t, args)
% Computes the normalized sinc function with bandwidth B at the point x

% t - time
% args = [B]
% B - bandwidth

if isa(args, 'struct')
    B = args.B;    
else
    B = args(1);  % For backwards compatibility
end

if t == 0,
    h = 2*B;
else
    h = sin(2*pi*B*t) / (pi*t);
end

% TESTING: delta function
%{
if abs(t) <= 0.001
    h = 5;
else
    h = 0;
end
%}

end
