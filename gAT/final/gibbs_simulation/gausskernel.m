function y = gausskernel(t,sigmah)

y=exp(-1/(2*sigmah^2).*t.^2);