function rmse = RMSE(v1, v2)
% Returns root mean squared error

rmse = sqrt(mean((v1 - v2).^2));
