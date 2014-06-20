function h = periodic_sinc(t, args)
% Periodic sinc function or Dirichlet kernel
% For Blu "Spares Sampling of Signal Innovations"

% t - time
% args has B and tau
% B - bandwidth
% tau - period

B = args.B;
tau = args.tau;

h = sin(pi*B*t)./(B*tau*sin(pi*t/tau));

%num = sin(pi*B*t)
%den = (B*tau*sin(pi*t/tau))

for i = 1 : length(t)
    if abs(sin(pi*t(i)/tau)) < 0.0001
        h(i) = 1;
    end
end

end
