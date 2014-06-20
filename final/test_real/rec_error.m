% Computes reconstruction error given
function e = rec_error(y, y_hat)

e = (norm(y_hat-y)/norm(y))^2;

end