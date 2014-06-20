% sinc kernel for test_cadzow_sinc
% t may be a vector

function h = sinc_kernel(t, B, tau)

h = sin(pi*B*t)./(B*tau*sin(pi*t/tau));

for i = 1 : length(t)
    if abs(sin(pi*t(i)/tau)) < 0.0001
        h(i) = 1;
    end
end

end