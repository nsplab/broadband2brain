% Script for testing SNR of channels
% Plots spike amplitude against noise RMSE
function SNR_vec = SNR_test(data_set, channels)

data_set = 1;
if nargin < 2
    channels = channels_to_use(data_set);
end

amp_vec = zeros(1, length(channels));
rmse_vec = zeros(1, length(channels));
SNR_vec = zeros(1, length(channels));

duration = 100;

for i = 1 : length(channels)
    
    % Data segment
    [training_start, training_end, data_start, data_end] = data_division(data_set);
    start_time = data_start;
    end_time = start_time + duration;

    % Get data
    buff = 0.1;
    [time, dat] = get_data(data_set, channels(i), start_time-buff, end_time+buff, 'raw');

    % Highpass
    [A, B] = butter_bp(300, 6000, 2);  % Change filter order?
    args.A_filter = A;
    args.B_filter = B;
    hp_handle = @filter_generic;  % Use causal version
    dat = -hp_handle(dat, args);
    
    rmse_vec(i) = sqrt(mean(dat.^2));
    amp_vec(i) = prctile(dat, 99.99);
    
    % Get real spikes
    real_spk = real_spikes(data_set, channels(i), start_time, end_time, 1);
    
    % special cases
    if channels(i) == 11
        amp_vec(i) = 50;
    end
    
    figure;
    hold on;
    num = 100000;
    plot(time(1:num), dat(1:num), 'c');
    xs = xlim;
    plot(xs, [rmse_vec(i) rmse_vec(i)], 'k--');
    plot(xs, [amp_vec(i) amp_vec(i)], 'k--');
    plot(real_spk, -10*ones(size(real_spk)), 'r.');
    xlim(xs)
    title(['channel ' int2str(channels(i)) ', ratio = ' num2str(amp_vec(i)/rmse_vec(i))]);
    
end

col = unifrnd(0, 1, length(channels), 3);

figure;
hold on;
for i = 1 : length(channels)
    plot(rmse_vec(i), amp_vec(i), 'LineStyle', 'none', 'Marker', '.', 'Color', col(i, :), 'MarkerSize', 20);
end
xlabel('RMSE');
ylabel('spike amplitude');
legend('3', '8', '9', '10', '11', '12', '13', '14', '15', '16');

SNR_vec = (amp_vec./rmse_vec).^2;

end