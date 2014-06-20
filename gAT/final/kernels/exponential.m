function h = exponential(t, args)
% Exponential decay, for calcium imaging
% Time constant tau

%if ~isfield(args, 'a')
%    args.a = 1
%end



%if t < 0,
%    h = 0;
%else
%    %h = args.a*exp(-t / args.tau);
%    h = exp(-t / args.tau);
%end

h = zeros(size(t));
h(t >= 0) = exp(-t(t >= 0) / args.tau);

end
