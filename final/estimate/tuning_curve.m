function [theta, freq, freq_b, pref_dir] = tuning_curve(vx, vy, spk, b)
% Returns theta, freq so that plot(theta, freq) gives tuning curve
% freq is (num times in theta bin and spike) / (num times in theta bin)
% freq_b is lambda based on b params
% pref_dir is preferred direction
% copied from preferred_dir.m

num_theta_bins = 20;

theta_spike = zeros(1, num_theta_bins);    % Number of times theta was in bin and spike
theta_total = zeros(1, num_theta_bins);    % Number of times theta was in bin

for i = 1 : length(vx),
    theta = atan2(vy(i), vx(i));
    theta_bin = ceil((theta + pi) * num_theta_bins / (2*pi));
    theta_total(theta_bin) = theta_total(theta_bin) + 1;
    if spk(i) >= 1,
        theta_spike(theta_bin) = theta_spike(theta_bin) + 1;
    end
end

theta_bins = ((1 : num_theta_bins) - 1/2) * 2*pi / num_theta_bins - pi;  % Centers of bins
preferred = theta_spike ./ theta_total;
vel = mean(sqrt(vx.^2 + vy.^2));
lambda = exp(b(1) + b(2) * vel * cos(theta_bins) + b(3) * vel * sin(theta_bins));

% Output
theta = theta_bins;
freq = preferred;
freq_b = lambda;
pref_dir = atan2(b(3), b(2));
pref_dir

end
