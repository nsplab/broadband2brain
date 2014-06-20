% Compare runtime of fast_oopsi and IterML

lambda = 1;

T_vec = [10 25 50 100 400];
dec_vec = [];
ml_vec = [];
for T = T_vec
    [runtime_dec runtime_ml] = test_cal(T, lambda);
    dec_vec = [dec_vec runtime_dec];
    ml_vec = [ml_vec runtime_ml];
end

% Plot
close all;
figure;
hold on;
plot(T_vec, log(dec_vec)./T_vec, 'r');
plot(T_vec, log(ml_vec)./T_vec, 'b');
legend('fast-oopsi (deconvolution)', 'IterML');
xlabel('T');
ylabel('runtime/T');