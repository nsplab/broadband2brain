% LFP summary plots
% scatter plots of LFP error metrics vs SNR
% plots of av LFP error metrics vs T (and grouped by SNR)
function[] = LFP_SNR()

close all;

dir = 'C:/Users/alex/Desktop/paper2plots/final/fig4/';

data_set = 1;
channels = channels_to_use(data_set);
SNR_vec = SNR_test(data_set, channels);
T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % Sampling periods (from uni_params)
T_ind = 6;

% scatter [absolute error in preferred phase] vs SNR for some specific T
T = T_vec(T_ind);
err_AT_lo = zeros(size(SNR_vec));
err_AT_hi = zeros(size(SNR_vec));
err_aFRI_lo = zeros(size(SNR_vec));
err_aFRI_hi = zeros(size(SNR_vec));
for i = 1 : length(channels)
    elec = channels(i);
    load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'pdlo_real', 'pdhi_real', 'pdlo_AT', 'pdhi_AT', 'pdlo_aFRI', 'pdhi_aFRI');
    err_AT_lo(i) = abs_err(pdlo_real, pdlo_AT(T_ind));
    err_AT_hi(i) = abs_err(pdhi_real, pdhi_AT(T_ind));
    err_aFRI_lo(i) = abs_err(pdlo_real, pdlo_aFRI(T_ind));
    err_aFRI_hi(i) = abs_err(pdhi_real, pdhi_aFRI(T_ind));
end
fig = figure;
hold on;
p1 = plot_dots(SNR_vec, err_AT_lo, lighten([1 0 0],1));
p2 = plot_dots(SNR_vec, err_AT_hi, [1 0 0]);
p3 = plot_dots(SNR_vec, err_aFRI_lo, lighten([0 0 1],1));
p4 = plot_dots(SNR_vec, err_aFRI_hi, [0 0 1]);
legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
title(['T = ' num2str(T)]);
xlabel('SNR');
ylabel('absolute error in preferred phase (degrees)');
print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\pref_SNR'], '-dpdf', '-r0');

% scatter KS p-value vs SNR for some specific T
T = T_vec(T_ind);
err_AT_lo = zeros(size(SNR_vec));
err_AT_hi = zeros(size(SNR_vec));
err_aFRI_lo = zeros(size(SNR_vec));
err_aFRI_hi = zeros(size(SNR_vec));
for i = 1 : length(channels)
    elec = channels(i);
    load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'KSlo_AT', 'KShi_AT', 'KSlo_aFRI', 'KShi_aFRI');
    err_AT_lo(i) = KSlo_AT(T_ind);
    err_AT_hi(i) = KShi_AT(T_ind);
    err_aFRI_lo(i) = KSlo_aFRI(T_ind);
    err_aFRI_hi(i) = KShi_aFRI(T_ind);
end
fig = figure;
hold on;
p1 = plot_dots(SNR_vec, err_AT_lo, lighten([1 0 0],1));
p2 = plot_dots(SNR_vec, err_AT_hi, [1 0 0]);
p3 = plot_dots(SNR_vec, err_aFRI_lo, lighten([0 0 1],1));
p4 = plot_dots(SNR_vec, err_aFRI_hi, [0 0 1]);
legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
title(['T = ' num2str(T)]);
xlabel('SNR');
ylabel('p-value (KS test)');
print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\KS_SNR'], '-dpdf', '-r0');

% scatter KS p-value vs SNR for some specific T
T = T_vec(T_ind);
err_AT_lo = zeros(size(SNR_vec));
err_AT_hi = zeros(size(SNR_vec));
err_aFRI_lo = zeros(size(SNR_vec));
err_aFRI_hi = zeros(size(SNR_vec));
for i = 1 : length(channels)
    elec = channels(i);
    load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'KSlo_AT', 'KShi_AT', 'KSlo_aFRI', 'KShi_aFRI');
    err_AT_lo(i) = KSlo_AT(T_ind);
    err_AT_hi(i) = KShi_AT(T_ind);
    err_aFRI_lo(i) = KSlo_aFRI(T_ind);
    err_aFRI_hi(i) = KShi_aFRI(T_ind);
end
fig = figure;
hold on;
p1 = plot_dots(SNR_vec, err_AT_lo, lighten([1 0 0],1));
p2 = plot_dots(SNR_vec, err_AT_hi, [1 0 0]);
p3 = plot_dots(SNR_vec, err_aFRI_lo, lighten([0 0 1],1));
p4 = plot_dots(SNR_vec, err_aFRI_hi, [0 0 1]);
legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
title(['T = ' num2str(T)]);
xlabel('SNR');
ylabel('p-value (KS test)');
print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\KS_SNR'], '-dpdf', '-r0');

% average[absolute error in preferred phase] vs T
runs = {'all channels', 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i)
        load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'pdlo_real', 'pdhi_real', 'pdlo_AT', 'pdhi_AT', 'pdlo_aFRI', 'pdhi_aFRI');
        err_AT_lo(i,:) = abs_err(pdlo_real, pdlo_AT);
        err_AT_hi(i,:) = abs_err(pdhi_real, pdhi_AT)
        err_aFRI_lo(i,:) = abs_err(pdlo_real, pdlo_aFRI);
        err_aFRI_hi(i,:) = abs_err(pdhi_real, pdhi_aFRI);
    end
    for i = 1 : length(T_vec)
        ci_AT_lo(:,i) = bootci(1000, @mean, err_AT_lo(:,i));
        ci_AT_hi(:,i) = bootci(1000, @mean, err_AT_hi(:,i));
        ci_aFRI_lo(:,i) = bootci(1000, @mean, err_aFRI_lo(:,i));
        ci_aFRI_hi(:,i) = bootci(1000, @mean, err_aFRI_hi(:,i));
    end
    fig = figure;
    hold on;
    %p1 = plot_errorbar(T_vec, mean(err_AT_lo), ci_AT_lo(1,:), ci_AT_lo(2,:), lighten([1 0 0],1));
    p2 = plot_errorbar(T_vec, mean(err_AT_hi), ci_AT_hi(1,:), ci_AT_hi(2,:), [1 0 0]);
    %p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), ci_aFRI_lo(1,:), ci_aFRI_lo(2,:), lighten([0 0 1],1));
    p4 = plot_errorbar(T_vec, mean(err_aFRI_hi), ci_aFRI_hi(1,:), ci_aFRI_hi(2,:), [0 0 1]);
    leg = legend([p2 p4], 'AT', 'gAT');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Sampling Period (T)', 'FontSize', 28);
    ylabel('Preferred Phase Error (degrees)', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim([0 0.053]);
    %title(runs{j});
    print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\pref_' names{j}], '-dpdf', '-r0');
    saveas(fig, [dir 'pref_' names{j}], 'epsc');
end

% average KS p-value vs T
runs = {'all channels', 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i);
        load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'KSlo_AT', 'KShi_AT', 'KSlo_aFRI', 'KShi_aFRI');
        err_AT_lo(i,:) = KSlo_AT;
        err_AT_hi(i,:) = KShi_AT;
        err_aFRI_lo(i,:) = KSlo_aFRI;
        err_aFRI_hi(i,:) = KShi_aFRI;
    end
    for i = 1 : length(T_vec)
        ci_AT_lo(:,i) = bootci(1000, @mean, err_AT_lo(:,i));
        ci_AT_hi(:,i) = bootci(1000, @mean, err_AT_hi(:,i));
        ci_aFRI_lo(:,i) = bootci(1000, @mean, err_aFRI_lo(:,i));
        ci_aFRI_hi(:,i) = bootci(1000, @mean, err_aFRI_hi(:,i));
    end
    fig = figure;
    hold on;
    p1 = plot_errorbar(T_vec, mean(err_AT_lo), ci_AT_lo(1,:), ci_AT_lo(2,:), lighten([1 0 0],1));
    p2 = plot_errorbar(T_vec, mean(err_AT_hi), ci_AT_hi(1,:), ci_AT_hi(2,:), [1 0 0]);
    p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), ci_aFRI_lo(1,:), ci_aFRI_lo(2,:), lighten([0 0 1],1));
    p4 = plot_errorbar(T_vec, mean(err_aFRI_hi), ci_aFRI_hi(1,:), ci_aFRI_hi(2,:), [0 0 1]);
    %leg = legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
    %set(leg, 'Location', 'NorthEast', 'FontSize', 24);
    xlabel('Sampling Period (T)', 'FontSize', 28);
    ylabel('p-Value (KS Test)', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim([0 0.053]);
    %title(runs{j});
    print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\KS_' names{j}], '-dpdf', '-r0');
    saveas(fig, [dir 'KS_' names{j}], 'epsc');
end

% average rayleigh p-value vs T
runs = {'all channels', 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i);
        load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'plo_AT', 'phi_AT', 'plo_aFRI', 'phi_aFRI');
        err_AT_lo(i,:) = plo_AT;
        err_AT_hi(i,:) = phi_AT;
        err_aFRI_lo(i,:) = plo_aFRI;
        err_aFRI_hi(i,:) = phi_aFRI;
    end
    for i = 1 : length(T_vec)
        ci_AT_lo(:,i) = bootci(1000, @mean, err_AT_lo(:,i));
        ci_AT_hi(:,i) = bootci(1000, @mean, err_AT_hi(:,i));
        ci_aFRI_lo(:,i) = bootci(1000, @mean, err_aFRI_lo(:,i));
        ci_aFRI_hi(:,i) = bootci(1000, @mean, err_aFRI_hi(:,i));
    end
    fig = figure;
    hold on;
    p1 = plot_errorbar(T_vec, mean(err_AT_lo), ci_AT_lo(1,:), ci_AT_lo(2,:), lighten([1 0 0],1));
    p2 = plot_errorbar(T_vec, mean(err_AT_hi), ci_AT_hi(1,:), ci_AT_hi(2,:), [1 0 0]);
    p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), ci_aFRI_lo(1,:), ci_aFRI_lo(2,:), lighten([0 0 1],1));
    p4 = plot_errorbar(T_vec, mean(err_aFRI_hi), ci_aFRI_hi(1,:), ci_aFRI_hi(2,:), [0 0 1]);
    %leg = legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
    %set(leg, 'Location', 'NorthEast', 'FontSize', 24);
    xlabel('Sampling Period (T)', 'FontSize', 28);
    ylabel('p-Value (Rayleigh Test)', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim([0 0.053]);
    %title(runs{j});
    print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\rayleigh_' names{j}], '-dpdf', '-r0');
    saveas(fig, [dir 'rayleigh_' names{j}], 'epsc');
end

% num rayleigh successes (alpha = 0.05) vs T
runs = {'all channels', 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
alpha = 0.05;
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i);
        load(['C:\Users\alex\Desktop\LFP\LFP_' int2str(data_set) '_' int2str(elec)], 'plo_AT', 'phi_AT', 'plo_aFRI', 'phi_aFRI');
        err_AT_lo(i,:) = plo_AT <= alpha;
        err_AT_hi(i,:) = phi_AT <= alpha;
        err_aFRI_lo(i,:) = plo_aFRI <= alpha;
        err_aFRI_hi(i,:) = phi_aFRI <= alpha;
    end
    fig = figure;
    hold on;
    p1 = plot_errorbar(T_vec, mean(err_AT_lo), 0, 0, lighten([1 0 0],1), 1);
    p2 = plot_errorbar(T_vec, mean(err_AT_hi), 0, 0, [1 0 0], 1);
    p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), 0, 0, lighten([0 0 1],1), 1);
    p4 = plot_errorbar(T_vec, mean(err_aFRI_hi), 0, 0, [0 0 1], 1);
    %leg = legend([p1 p2 p3 p4], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)');
    %set(leg, 'Location', 'NorthEast', 'FontSize', 24);
    xlabel('Sampling Period (T)', 'FontSize', 28);
    ylabel('Proportion of Rayleigh Successes', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim([0 0.053]);
    %title(runs{j});
    print(fig, ['C:\Users\alex\Desktop\paper2plots\LFP\summary\rayleigh_decision_' names{j}], '-dpdf', '-r0');
    saveas(fig, [dir 'rayleigh_frac_' names{j}], 'epsc');
end

end

% compute absolute error in direction (degrees)
% if other is a vector, computes vector of results
function err = abs_err(real, other)
err = abs(real-other);
err = min(err, 360-err);
end