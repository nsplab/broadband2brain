% Plot true positive rate vs false positive rate
% Params specified in test_generic_reconstruct

mult_vec = [-3 : 0.1 : 3];
tp_vec = zeros(size(mult_vec));
fp_vec = zeros(size(mult_vec));

paramID = 34;

for i = 1 : length(mult_vec)
    [samples, samples_model, tp_rate, fp_rate] = test_generic_reconstruct(paramID, mult_vec(i));
    tp_vec(i) = tp_rate;
    fp_vec(i) = fp_rate;
end

close all;
figure;
plot(tp_vec, fp_vec);
xlabel('true positive rate');
ylabel('false positive rate');

tp_vec
fp_vec