% Plot false negatives vs false positives

close all;

methods = {'ann', 'gibbs'};%, 'RMSE'};

data_set = 1;
electrodes = 9;%channels_to_use(1);

colors = {[0 0 0.5], [0 0.5 0], [1 0 1]};

for electrode = electrodes

    figure;
    hold on;
    
    h = [];

    for i = 1 : length(methods)
        
        % 400 Hz
        %sr = 400;
        paramIDs = [11, 51, 11];
        
        [fn_vec fp_vec] = tp_fp_2(methods{i}, paramIDs(i), data_set, electrode);
        h1 = plot_errorbar(fn_vec, fp_vec, fp_vec, fp_vec, colors{i});
        h = [h h1];

        % 1000 Hz
        %sr = 1000;
        paramIDs = [14, 54, 14];
        
        [fn_vec fp_vec] = tp_fp_2(methods{i}, paramIDs(i), data_set, electrode);
        h1 = plot_errorbar(fn_vec, fp_vec, fp_vec, fp_vec, lighten(colors{i}));
        h = [h h1];

        % 3000 Hz
        %sr = 3000
        paramIDs = [17, 57, 17];
        
        [fn_vec fp_vec] = tp_fp_2(methods{i}, paramIDs(i), data_set, electrode);
        h1 = plot_errorbar(fn_vec, fp_vec, fp_vec, fp_vec, lighten(lighten(colors{i})));
        h = [h h1];
    end

    legend(h, 'annihilating filter (400 Hz)', 'annihilating filter (1000 Hz)', 'annihilating filter (3000 Hz)', ...
        'gibbs (400 Hz)', 'gibbs (1000 Hz)', 'gibbs (3000 Hz)');

    xlabel('missed spikes (as percent of total true spikes)');
    ylabel('false positives per second');
    title(['false positives vs missed spikes, channel ' int2str(electrode)]);
    
    ylim([0 100]);
    
end

saveas(f, ['C:\Users\alex\Desktop\fn_fp'], 'epsc');