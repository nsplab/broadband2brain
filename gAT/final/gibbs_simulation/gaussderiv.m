% Derivative of Gaussian Kernel
% See gausskernel.m
% Used for computation of Cramer Rao bound

function y = gaussderiv(t, sigmah)

y = -t/sigmah^2.*exp(-1/(2*sigmah^2).*t.^2);

end