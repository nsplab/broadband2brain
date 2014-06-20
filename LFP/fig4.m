% Fig 4: LFP
function fig4()

close all;

dir = [root_dir() '../paper2plots/final/fig4/'];

data_set = 1;
channels = channels_to_use(data_set);
T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods (from uni_params)

% average[absolute error in preferred phase] vs T
runs = {'all channels'};%, 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    err_TD_lo = zeros(length(chan), length(T_vec));
    err_TD_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i);
        load([root_dir() '../LFP/LFP_' int2str(data_set) '_' int2str(elec)], 'pdlo_real', 'pdhi_real', 'pdlo_AT', 'pdhi_AT', 'pdlo_aFRI', 'pdhi_aFRI', 'pdlo_TD', 'pdhi_TD', 'pdglm_real', 'pdglm_AT', 'pdglm_aFRI', 'pdglm_TD');
        %err_AT_lo(i,:) = abs_err(pdlo_real, pdlo_AT);
        err_AT_hi(i,:) = abs_err(pdglm_real, pdglm_AT);
        %err_aFRI_lo(i,:) = abs_err(pdlo_real, pdlo_aFRI);
        err_aFRI_hi(i,:) = abs_err(pdglm_real, pdglm_aFRI);
        %err_TD_lo(i,:) = abs_err(pdlo_real, pdlo_TD);
        err_TD_hi(i,:) = abs_err(pdglm_real, pdglm_TD);
    end
    for i = 1 : length(T_vec)
        %ci_AT_lo(:,i) = bootci(1000, @mean, err_AT_lo(:,i));
        ci_AT_hi(:,i) = bootci(1000, @mean, err_AT_hi(:,i));
        %ci_aFRI_lo(:,i) = bootci(1000, @mean, err_aFRI_lo(:,i));
        ci_aFRI_hi(:,i) = bootci(1000, @mean, err_aFRI_hi(:,i));
        %ci_TD_lo(:,i) = bootci(1000, @mean, err_TD_lo(:,i));
        ci_TD_hi(:,i) = bootci(1000, @mean, err_TD_hi(:,i));
    end
    fig = figure;
    hold on;
    %p1 = plot_errorbar(T_vec, mean(err_AT_lo), ci_AT_lo(1,:), ci_AT_lo(2,:), lighten([1 0 0],1));
    p2 = plot_errorbar(1000*T_vec, mean(err_AT_hi), ci_AT_hi(1,:), ci_AT_hi(2,:), [1 0 0]);
    %p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), ci_aFRI_lo(1,:), ci_aFRI_lo(2,:), lighten([0 0 1],1));
    p4 = plot_errorbar(1000*T_vec, mean(err_aFRI_hi), ci_aFRI_hi(1,:), ci_aFRI_hi(2,:), [0 0 1]);
    %p5 = plot_errorbar(T_vec, mean(err_TD_lo), ci_TD_lo(1,:), ci_TD_lo(2,:), lighten([0 1 0],1));
    p6 = plot_errorbar(1000*T_vec, mean(err_TD_hi), ci_TD_hi(1,:), ci_TD_hi(2,:), [0 1 0]);
    leg = legend([p2 p4 p6], 'AT', 'gAT', 'twodelta');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Preferred Phase Error (degrees)', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim(1000*[0 0.103]);
    %title(runs{j});
    %print(fig, 'C:\Users\alex\Desktop\paper2plots\final\fig4\fig4b', '-dpdf', '-r0');
    saveas(fig, [dir 'fig4b'], 'epsc');
    saveas(fig, [dir 'fig4b'], 'fig');
end

% av absolute error in distribution vs T
runs = {'all channels'};%, 'low SNR', 'high SNR'};
names = {'all', 'low', 'high'};
chans = {1:10, [4 5 6 7 8], [1 2 3 9 10]};  % indices in channels
for j = 1 : length(runs)
    chan = channels(chans{j});  % which channels to average over
    err_AT_lo = zeros(length(chan), length(T_vec));
    err_AT_hi = zeros(length(chan), length(T_vec));
    err_aFRI_lo = zeros(length(chan), length(T_vec));
    err_aFRI_hi = zeros(length(chan), length(T_vec));
    err_TD_lo = zeros(length(chan), length(T_vec));
    err_TD_hi = zeros(length(chan), length(T_vec));
    for i = 1 : length(chan)
        elec = chan(i);
        load([root_dir() '../LFP/LFP_' int2str(data_set) '_' int2str(elec)], 'distrhi_real', 'distrhi_AT', 'distrhi_aFRI', 'distrhi_TD');
        %err_AT_lo(i,:) = abs_err(pdlo_real, pdlo_AT);
        err_AT_hi(i,:) = mean(abs(repmat(distrhi_real',[length(T_vec) 1])-distrhi_AT)');
        %err_aFRI_lo(i,:) = abs_err(pdlo_real, pdlo_aFRI);
        err_aFRI_hi(i,:) = mean(abs(repmat(distrhi_real',[length(T_vec) 1])-distrhi_aFRI)');
        %err_TD_lo(i,:) = abs_err(pdlo_real, pdlo_TD);
        err_TD_hi(i,:) = mean(abs(repmat(distrhi_real',[length(T_vec) 1])-distrhi_TD)');
    end
    for i = 1 : length(T_vec)
        %ci_AT_lo(:,i) = bootci(1000, @mean, err_AT_lo(:,i));
        ci_AT_hi(:,i) = bootci(1000, @mean, err_AT_hi(:,i));
        %ci_aFRI_lo(:,i) = bootci(1000, @mean, err_aFRI_lo(:,i));
        ci_aFRI_hi(:,i) = bootci(1000, @mean, err_aFRI_hi(:,i));
        %ci_TD_lo(:,i) = bootci(1000, @mean, err_TD_lo(:,i));
        ci_TD_hi(:,i) = bootci(1000, @mean, err_TD_hi(:,i));
    end
    fig = figure;
    hold on;
    %p1 = plot_errorbar(T_vec, mean(err_AT_lo), ci_AT_lo(1,:), ci_AT_lo(2,:), lighten([1 0 0],1));
    p2 = plot_errorbar(1000*T_vec, mean(err_AT_hi), ci_AT_hi(1,:), ci_AT_hi(2,:), [1 0 0]);
    %p3 = plot_errorbar(T_vec, mean(err_aFRI_lo), ci_aFRI_lo(1,:), ci_aFRI_lo(2,:), lighten([0 0 1],1));
    p4 = plot_errorbar(1000*T_vec, mean(err_aFRI_hi), ci_aFRI_hi(1,:), ci_aFRI_hi(2,:), [0 0 1]);
    %p5 = plot_errorbar(T_vec, mean(err_TD_lo), ci_TD_lo(1,:), ci_TD_lo(2,:), lighten([0 1 0],1));
    p6 = plot_errorbar(1000*T_vec, mean(err_TD_hi), ci_TD_hi(1,:), ci_TD_hi(2,:), [0 1 0]);
    leg = legend([p2 p4 p6], 'AT', 'gAT', 'twodelta');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Av. Distribution Error', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim(1000*[0 0.103]);
    %title(runs{j});
    %print(fig, 'C:\Users\alex\Desktop\paper2plots\final\fig4\fig4c', '-dpdf', '-r0');
    saveas(fig, [dir 'fig4c'], 'epsc');
    saveas(fig, [dir 'fig4c'], 'fig');
end

% baseline firing rate scatter plot
% gAT baseline vs AT baseline vs twodelta baseline
%{
% TODO: ??? HOW TO DO THIS...
chans = 1:10;
T_scatter = 0.005;  % can change this
methods = {'real', 'AT', 'gAT', 'twodelta'};
for j = 1 : length(runs)
    %vec_real = zeros(1, length(chans));
    vec_AT = zeros(1, length(chans));
    vec_gAT = zeros(1, length(chans));
    for i = 1 : length(chans)
        elec = chan(i);
        % AT
        method = 'AT';
        load([root_dir() '../LFP/distr/distr_' method '_' int2str(data_set) '_' int2str(elec) '_' int2str(1000*T_scatter)]);
        vx = cos(lfp_vec); % not actually velocities but analogous to velocities in tuning curve calculation
        vy = sin(lfp_vec);
        [b, dev, stats] = glmfit([vx vy], spk_vec, 'poisson', 'const', 'on');
        vec_AT(i) = b(1);
        % gAT
        method = 'gAT';
        load([root_dir() '../LFP/distr/distr_' method '_' int2str(data_set) '_' int2str(elec) '_' int2str(1000*T_scatter)]);
        vx = cos(lfp_vec); % not actually velocities but analogous to velocities in tuning curve calculation
        vy = sin(lfp_vec);
        [b, dev, stats] = glmfit([vx vy], spk_vec, 'poisson', 'const', 'on');
        vec_gAT(i) = b(1);
    end
    fig = figure;
    hold on;
    scatter(vec_AT, vec_gAT, 'ro');
    plot([-8 -2], [-8 -2], 'k');
    xlabel('Baseline Firing Rate, AT', 'FontSize', 28);
    ylabel('Baseline Firing Rate, gAT', 'FontSize', 28);
    set(gca, 'FontSize', 24);
    xlim([-8 -2]);
    ylim([-8 -2]);
    %title(runs{j});
    %print(fig, [root_dir() '../paper2plots/final/fig4/fig4c', '-dpdf', '-r0');
    saveas(fig, [dir 'fig4_scatter'], 'epsc');
end
%}

end

% compute absolute error in direction (degrees)
% if other is a vector, computes vector of results
function err = abs_err(real, other)
err = abs(real-other);
err = min(err, 360-err);
end
