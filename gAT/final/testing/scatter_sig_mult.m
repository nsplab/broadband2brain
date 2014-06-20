% Looks for best value of multiple of sigma for rejection threshold

method_name = 'int';
paramID = 1;

std_dev = [];
mult = [];

for data_set = [1, 2];

    for elec = channels_to_use(data_set)

        [sig, best_multiple] = test_stddev_thres(method_name, paramID, data_set, elec);
        std_dev = [std_dev sig];
        mult = [mult best_multiple];

    end

end

close all;

figure;
plot(std_dev, mult, 'b.');
xlabel('std dev');
ylabel('multiple');

figure;
plot(std_dev, mult .* std_dev, 'b.');
xlabel('std dev');
ylabel('thres');

mult_ind = find(mult < 5);
mult = mult(mult_ind);
av_mult = mean(mult);
