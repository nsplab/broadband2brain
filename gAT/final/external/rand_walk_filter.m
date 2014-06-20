function [px_r, py_r, vx_r, vy_r, ax_r, ay_r] = rand_walk_filter(K, C, dt, start, W_0, A, Q, spike_vec, b_0, b_1, b_2)

% Point process filter with random walk prior
% Model: theta_{k+1} = A * theta_k + w_k where w_k has covariance matrix Q

% Author: Alex Wein, July 2011

% Inputs:
% K - number of timesteps
% C - number of neurons
% start - start state (column vector)
% W_0 - variance on start state
% A - state update matrix
% Q - random walk covariance
% spike_vec - spike_vec(c, i) = 0 or 1 depending on if neuron c has spike at timestep i
% b_0, b_1, b_2 vectors of length C containing model parameters

% Outputs:
% Decoded px, py, vx, vy

num_dims = size(W_0, 1);  % Length of state vector

theta = zeros(K+1, num_dims);  % Vector of px, py, vx, vy
W = zeros(K+1, num_dims, num_dims);   % Variance

% Initial conditions k = 0
theta(1, :) = start';
W(1, :, :) = W_0;

for k = 1 : K,
    %disp(['k=' int2str(k) '  dN_total=' int2str(sum(spike_vec(:, k)))]);
    theta(k+1, :) = A * theta(k, :)';  % (A.1)
    W(k+1, :, :) = A * squeeze(W(k, :, :)) * A' + Q;  % (A.2)

    sum1 = zeros(num_dims, 1);  % for (A.7)
    sum2 = zeros(num_dims);  % for (A.8)
    for c = 1 : C,
        dN = spike_vec(c, k);
        lambda = exp(b_0(c) + b_1(c) * theta(k+1, 3) + b_2(c) * theta(k+1, 4));
        %disp(['max lambda (per dt): ' num2str(max(lambda) * dt)]);
        vec = zeros(num_dims, 1);
        vec(3) = b_1(c);
        vec(4) = b_2(c);
        sum1 = sum1 + vec * (dN - lambda*dt);
        sum2 = sum2 + vec * vec' * lambda * dt;
    end

    W(k+1, :, :) = (squeeze(W(k+1, :, :))^(-1) + sum2)^(-1);  % (A.8)
    theta(k+1, :) = theta(k+1, :)' + squeeze(W(k+1, :, :)) * sum1;  % (A.7)
end

px_r = theta(:, 1);
py_r = theta(:, 2);
vx_r = theta(:, 3);
vy_r = theta(:, 4);
ax_r = 0;
ay_r = 0;
if size(theta, 2) >= 6,
    ax_r = theta(:, 5);
    ay_r = theta(:, 6);
end
