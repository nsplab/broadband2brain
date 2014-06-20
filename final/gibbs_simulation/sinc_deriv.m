% Derivative of periodic sinc kernel
% t may be a vector

function h = sinc_deriv(t, B, tau)

h = (pi*B*tau*cos(pi*B*t)*sin(pi*t/tau)-pi*sin(pi*B*t)*cos(pi*t/tau))/(B*tau^2*sin(pi*t/tau)^2);

for i = 1 : length(t)
    if abs(sin(pi*t(i)/tau)) < 0.0001
        h = 0;
    end
end

end