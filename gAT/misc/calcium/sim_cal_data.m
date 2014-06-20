% Simulate calcium imaging data

% T = num time steps
% delta = length of time step (sec)
% lambda = av firing rate (Hz)
% alpha, beta, sigma

% F_t ~ N(alpha*C_t + beta, sigma^2)
% P[n_t] = Poisson(lambda*delta)
% C_t = gamma*C_{t-1} + n_t
% gamma = (1 - delta*tau)

function [n C F] = sim_cal_data(T, delta, lambda, alpha, beta, sigma, gamma)

n = poissrnd(lambda*delta, 1, T);
C = zeros(size(n));
F = zeros(size(n));

for t = 2 : length(C)
    C(t) = gamma*C(t-1) + n(t);
    F(t) = normrnd(alpha*C(t) + beta, sigma);
end

end