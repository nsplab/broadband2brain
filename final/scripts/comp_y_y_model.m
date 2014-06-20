% Compare samples with model-reconstructed samples
% As a function of sampling rate
% Hoping to explain why performance degrades at high sampling rates

close all;
clear all;

% copied from uni_params.m (removed 200 Hz)
sampling_rates = [400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
paramIDs = 11:19;

% add low sampling rates
sampling_rates = [20, 40, 60, 80, 100, 120, 150, 200, 250, 300 sampling_rates];
paramIDs = [31:40 paramIDs];

RMSE_vec = zeros(size(paramIDs));

for i = 1 : length(paramIDs)
    [y, y_model] = test_generic_reconstruct(paramIDs(i));
    RMSE_vec(i) = RMSE(y, y_model);
end
    
figure;
plot(sampling_rates, RMSE_vec);
xlabel('sampling rate (Hz)');
ylabel('RMSE');
title('RMSE between samples and model-reconstructed samples');