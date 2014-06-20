% Plot error in preferred direction vs sampling rate for various methods

%%{
sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000];%, 3000, 4000, 5000];
paramIDs = 0 : 6;

method_names = {'ann', 'gibbs'};%, 'RMSE'};
paramOffset = [10, 30, 10];
data_set = 1;
electrodes = channels_to_use(data_set);

colors = {[0 0 1], [0 1 0], [1 0 1]};

h = [];

data = [];

for i = 1 : length(method_names)
    vec = [];
    for j = 1 : length(sampling_rates)
        [real_dir_vec dir_vec dir_err] = est_delay_b(method_names{i}, paramIDs(j)+paramOffset(i), data_set, 1, electrodes);
        vec(j) = dir_err;
    end
    data = [data ; vec];
end

save('C:\Users\alex\Desktop\figures\mat\pref');
%%}
%load('C:\Users\alex\Desktop\figures\mat\pref');

close all;
f = figure;
hold on;

data = data * 180 / pi;

for i = 1 : 2%length(methods)
    h1 = plot_errorbar(sampling_rates, data(i,:), data(i,:), data(i,:), colors{i});
    h = [h h1];
end

legend(h, 'annihilating filter', 'gibbs');%, 'IterML');

xlabel('sampling rate (Hz)');
ylabel('average absolute error in preferred direction (degrees)');
title('error in preferred direction vs sampling rate, subject 1');

ylim([0 180]);

saveas(f, ['C:\Users\alex\Desktop\prefdir'], 'epsc');