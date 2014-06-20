% Plot t_error vs T (for each channel)
% Plot fraction intervals vs T (for each channel)
% Fraction of intervals: of intervals with >=1 real spikes, how many have
% exactly 1 real and 1 reconstructed
% Uses thresholds from get_thres.m

close all;
clear all;

load = 1;
load

if load == 0

    % Parameters
    data_set = 1;
    channels = channels_to_use(data_set);
    duration = 60;
    start_offset = 60;
    T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods
    
    te_AT = cell(16, length(T_vec));  % t_AT(elec, T_ind) = list of abs(t_error) values
    te_aFRI = cell(16, length(T_vec));
    te_TD = cell(16, length(T_vec));
    num_AT = zeros(16, length(T_vec));  % num_AT(elec, T_ind) = number of intervals used (numerator)
    num_aFRI = zeros(16, length(T_vec));
    num_TD = zeros(16, length(T_vec));
    num_TD1 = zeros(16, length(T_vec));
    num_TD2 = zeros(16, length(T_vec));
    den_AT = zeros(16, length(T_vec));  % den_aFRI(elec_T_ind) = number of intervals with spikes (denominator)
    den_aFRI = zeros(16, length(T_vec));
    den_TD = zeros(16, length(T_vec));
    den_TD1 = zeros(16, length(T_vec));
    den_TD2 = zeros(16, length(T_vec));
    
    for i = 1 : 16
        for j = 1 : length(T_vec)
            t_AT{i, j} = [];
            t_aFRI{i, j} = [];
            t_TD{i, j} = [];
        end
    end

    % Loop through methods
    methods = {'Analog Thresholding', 'gAT', 'twodelta'};
    colors = {[1 0 0], [0 0 1], [0 1 0]};


    % Count the number of trials per method
    numTrials = 0;
    for q = 1 : length(channels)
        % Data segment
        [training_start, training_end, data_start, data_end] = data_division(data_set);
        start_time = data_start;
        end_time = start_time + duration;

        % Loop through T
        for i = 1 : length(T_vec)
            T = T_vec(i);
            t_block = start_time : T : end_time;
            t_block = t_block(1 : end-1);  % t_block contains start times
            % Loop through intervals
            for k = 1 : length(t_block)
                numTrials = numTrials + 1;
            end
        end
    end

    % Initialize the lists for storing the output
    Tlist = zeros(numTrials, 1);
    ylist = zeros(numTrials, 4);
    yguess3list = zeros(numTrials, 1);
    yguess4list = zeros(numTrials, 1);
    num_rec_spikes_TDlist = zeros(numTrials, 1);
    num_real_spikeslist = zeros(numTrials, 1);
    trialCount = 0;
    num_rec_spikes_Flist = zeros(numTrials, 1);
    trialCount2 = 0;
    
    for q = 1 : length(channels)
        elec = channels(q);
        electrode = elec
        fprintf(sprintf('q: %d / %d\n', q, length(channels)));
    
        % Data segment
        [training_start, training_end, data_start, data_end] = data_division(data_set);
        start_time = data_start;
        end_time = start_time + duration;

        % Get data
        buff = 0.1;
        [time, dat] = get_data(data_set, electrode, start_time-buff, end_time+buff, 'raw');

        % Highpass
        [A, B] = butter_bp(300, 6000, 2);  % Change filter order?
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        dat = -hp_handle(dat, args);

        % Get real spikes
        real_spk = real_spikes(data_set, electrode, start_time, end_time, 1);
        %av_spk_rate = length(real_spk)/duration
        
        % Load thres
        thres = get_thres(data_set, elec);

        for m = 1 : length(methods)
            fprintf(sprintf('m: %d / %d\n', m, length(methods)));

            % Loop through T
            for i = 1 : length(T_vec)
                fprintf(sprintf('i: %d / %d\n', i, length(T_vec)));
                T = T_vec(i);
                t_block = start_time : T : end_time;
                t_block = t_block(1 : end-1);  % t_block contains start times

                real_ind = 1;  % index in real_spk

                t_rec = [];

                % Loop through intervals
                for k = 1 : length(t_block)

                    % find real spikes in interval
                    num_real_spikes = 0;
                    t_real = [];
                    while real_ind <= length(real_spk) && real_spk(real_ind) <= t_block(k)+T
                        t_real = [t_real real_spk(real_ind) - t_block(k)];
                        num_real_spikes = num_real_spikes + 1;
                        real_ind = real_ind + 1;
                    end

                    % get data segment
                    seg_ind = find_fast(time, t_block(k), t_block(k)+T, get_dt());
                    seg = dat(seg_ind);

                    if m == 1
                        % get reconstructed spikes in interval
                        % analog thresholding
                        t_AT = T/2;
                        num_rec_spikes_AT = max(seg) > thres;

                        t = t_AT;
                        num_rec_spikes = num_rec_spikes_AT;
                    elseif m == 2
                        % aFRI
                        seg_comp = seg > thres;
                        y = succ_int(seg_comp, get_dt(), 2);
                        [t_F, c] = reconstruct_analog(y, T, [], []);
                        num_rec_spikes_F = (c ~= 0);

                        t = t_F;
                        num_rec_spikes = num_rec_spikes_F;

                        trialCount2 = trialCount2 + 1;
                        num_rec_spikes_Flist(trialCount2, :) = num_rec_spikes_F;
                    elseif m == 3
                        % twodelta
                        seg_comp = seg > thres;
                        y = succ_int(seg_comp, get_dt(), 4);
                        [t_TD, c, sigma, elapsed_time, yguess3, yguess4] = reconstruct_twodelta(y, T, [], []);
                        num_rec_spikes_TD = sum(c ~= 0);

                        t = t_TD;
                        num_rec_spikes = num_rec_spikes_TD;

                        trialCount = trialCount + 1;
                        Tlist(trialCount, :) = T;
                        ylist(trialCount, :) = y';
                        yguess3list(trialCount, :) = yguess3;
                        yguess4list(trialCount, :) = yguess4;
                        num_rec_spikes_TDlist(trialCount, :) = num_rec_spikes_TD;
                        num_real_spikeslist(trialCount, :) = num_real_spikes;
                    end


                    if m == 1
                        if num_real_spikes > 0
                            den_AT(elec, i) = den_AT(elec, i) + 1;
                            if num_rec_spikes == num_real_spikes
                                num_AT(elec, i) = num_AT(elec, i) + 1;
                            end
                        end
                        if num_rec_spikes == 1 && num_real_spikes == 1
                            te_AT{elec, i} = [te_AT{elec, i} abs(t - t_real)];
                        end
                    elseif m == 2
                        if num_real_spikes > 0
                            den_aFRI(elec, i) = den_aFRI(elec, i) + 1;
                            if num_rec_spikes == num_real_spikes
                                num_aFRI(elec, i) = num_aFRI(elec, i) + 1;
                            end
                        end
                        if num_rec_spikes == 1 && num_real_spikes == 1
                            te_aFRI{elec, i} = [te_aFRI{elec, i} abs(t - t_real)];
                        end
                    elseif m == 3
                        if num_real_spikes > 0
                            den_TD(elec, i) = den_TD(elec, i) + 1;
                            den_TD1(elec, i) = den_TD1(elec, i) + (num_real_spikes == 1);
                            den_TD2(elec, i) = den_TD2(elec, i) + (num_real_spikes == 2);
                            if num_rec_spikes == num_real_spikes
                                num_TD(elec, i) = num_TD(elec, i) + 1;
                                num_TD1(elec, i) = num_TD1(elec, i) + (num_real_spikes == 1);
                                num_TD2(elec, i) = num_TD2(elec, i) + (num_real_spikes == 2);
                            end
                        end
                        if num_rec_spikes == 1 && num_real_spikes == 1
                            te_TD{elec, i} = [te_TD{elec, i} abs(t(1) - t_real)];
                        end
                    end

                end

            end

        end

    end

    save([root_dir() 'paper2/mat/fig2']);

else
    
    clear all;
    load([root_dir() 'paper2/mat/fig2']);
end
    
dir = [root_dir() '../paper2plots/final/fig2/'];

% Pick one of these:

channels_to_plot = channels;
name = 'all';
do_mean = 1;

%channels_to_plot = [14];
%name = 'single';
%do_mean = 0;

% t_error vs T (separate channels)
fig = figure;
hold on;
data_AT = zeros(1, length(T_vec));
data_aFRI = zeros(1, length(T_vec));
data_TD = zeros(1, length(T_vec));
cihi_AT = zeros(1, length(T_vec));
hilo_AT = zeros(1, length(T_vec));
cihi_aFRI = zeros(1, length(T_vec));
cilo_aFRI = zeros(1, length(T_vec));
cihi_TD = zeros(1, length(T_vec));
cilo_TD = zeros(1, length(T_vec));
for elec = channels_to_plot
    for i = 1 : length(T_vec)
        data_AT(i) = mean(te_AT{elec, i});
        ci = bootci(1000, @mean, te_AT{elec, i});
        cilo_AT(i) = ci(1);
        cihi_AT(i) = ci(2);
        data_aFRI(i) = mean(te_aFRI{elec, i});
        ci = bootci(1000, @mean, te_aFRI{elec, i});
        cilo_aFRI(i) = ci(1);
        cihi_aFRI(i) = ci(2);
        if (size(te_TD{elec, i}) ~= [0 0])
            data_TD(i) = mean(te_TD{elec, i});
            ci = bootci(1000, @mean, te_TD{elec, i});
            cilo_TD(i) = ci(1);
            cihi_TD(i) = ci(2);
        else
            data_TD(i) = NaN;
            cilo_TD(i) = NaN;
            cihi_TD(i) = NaN;
        end
    end
    p1 = plot_errorbar(1000*T_vec, 1000*data_AT, 0, 0, [1 0 0], 1);
    p2 = plot_errorbar(1000*T_vec, 1000*data_aFRI, 0, 0, [0 0 1], 1);
    p3 = plot_errorbar(1000*T_vec, 1000*data_TD, 0, 0, [0 1 0], 1);
end
xlim(1000*[0 1.06*max(T_vec)]);
ylim(1000*[0 2*0.014]);
p0 = plot(xlim, 0.25*xlim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Av. Spike Time Error (ms)', 'FontSize', 28);
leg = legend([p0 p1 p2 p3], 'AT theoretical', 'AT', 'gAT', 'twodelta');
set(leg, 'Location', 'NorthWest', 'FontSize', 28);
set(gca, 'FontSize', 24);
saveas(fig, [dir 't_error_' name], 'epsc');
saveas(fig, [dir 't_error_' name], 'fig');

% frac used vs T (separate channels)
fig = figure;
hold on;
data_AT = zeros(1, length(T_vec));
data_aFRI = zeros(1, length(T_vec));
data_TD = zeros(1, length(T_vec));
for elec = channels_to_plot
    for i = 1 : length(T_vec)
        data_AT(i) = num_AT(elec, i)/den_AT(elec,i);
        %data_aFRI(i) = num_aFRI(elec, i)/den_aFRI(elec,i);
        data_TD(i) = num_TD(elec, i)/den_TD(elec,i);
    end
    p1 = plot_errorbar(1000*T_vec, data_AT, 0, 0, [1 0 0], 1);
    %p2 = plot_errorbar(1000*T_vec, data_aFRI, 0, 0, [0 0 1], 1);
    p3 = plot_errorbar(1000*T_vec, data_TD, 0, 0, [0 1 0], 1);
end
leg = legend([p1 p3], 'gAT', 'twodelta', 'Location', 'Best');
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Fraction of Valid Intervals', 'FontSize', 28);
%leg = legend([p1], 'Both methods');
%set(leg, 'Location', 'NorthWest', 'FontSize', 28);
xlim(1000*[0 1.06*max(T_vec)]);
ylim([0 1]);
set(gca, 'FontSize', 24);
saveas(fig, [dir 'frac_' name], 'epsc');
saveas(fig, [dir 'frac_' name], 'fig');

% frac valid for TD (1 spike, 2 spike)
fig = figure;
hold on;
data_TD1 = zeros(1, length(T_vec));
data_TD2 = zeros(1, length(T_vec));
for elec = channels_to_plot
    for i = 1 : length(T_vec)
        data_TD1(i) = num_TD1(elec, i)/den_TD1(elec,i);
        data_TD2(i) = num_TD2(elec, i)/den_TD2(elec,i);
    end
    p1 = plot_errorbar(1000*T_vec, data_TD1, 0, 0, [1 0 0], 1);
    p2 = plot_errorbar(1000*T_vec(~isnan(data_TD2)), data_TD2(~isnan(data_TD2)), 0, 0, [0 1 0], 1);
    data_TD2
end
leg = legend([p1 p2], '1 spike', '2 spike', 'Location', 'Best');
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Fraction of Valid Intervals', 'FontSize', 28);
xlim(1000*[0 1.06*max(T_vec)]);
ylim([0 1]);
set(gca, 'FontSize', 24);
saveas(fig, [dir 'frac_TD12_' name], 'epsc');
saveas(fig, [dir 'frac_TD12_' name], 'fig');

% difference in frac used vs T (separate channels)
fig = figure;
hold on;
data_AT = zeros(1, length(T_vec));
data_aFRI = zeros(1, length(T_vec));
data_TD = zeros(1, length(T_vec));
for elec = channels_to_plot
    for i = 1 : length(T_vec)
        data_AT(i) = num_AT(elec, i)/den_AT(elec,i);
        %data_aFRI(i) = num_aFRI(elec, i)/den_aFRI(elec,i);
        data_TD(i) = num_TD(elec, i)/den_TD(elec,i);
    end
    p1 = plot_errorbar(1000*T_vec, data_TD - data_AT, 0, 0, [1 0 0], 1);
end
leg = legend([p1], 'Fraction Valid of TD - Fraction Valid of AT/gAT', 'Location', 'Best');
xlabel('Sampling Period T (ms)', 'FontSize', 28);
ylabel('Difference in Fraction of Valid Intervals', 'FontSize', 28);
%leg = legend([p1], 'Both methods');
%set(leg, 'Location', 'NorthWest', 'FontSize', 28);
xlim(1000*[0 1.06*max(T_vec)]);
%ylim([0 1]);
set(gca, 'FontSize', 24);
saveas(fig, [dir 'diff_frac_' name], 'epsc');
saveas(fig, [dir 'diff_frac_' name], 'fig');

if do_mean

    % t_error vs T (mean + bootci)
    fig = figure;
    hold on;
    data_AT = zeros(1, length(T_vec));
    data_aFRI = zeros(1, length(T_vec));
    data_TD = zeros(1, length(T_vec));
    cihi_AT = zeros(1, length(T_vec));
    hilo_AT = zeros(1, length(T_vec));
    cihi_aFRI = zeros(1, length(T_vec));
    cilo_aFRI = zeros(1, length(T_vec));
    cihi_TD = zeros(1, length(T_vec));
    cilo_TD = zeros(1, length(T_vec));
    for i = 1 : length(T_vec)
        vec = [];
        for elec = channels_to_plot
            vec = [vec mean(te_AT{elec, i})];
        end
        data_AT(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_AT(i) = ci(1);
        cihi_AT(i) = ci(2);
        vec = [];
        for elec = channels_to_plot
            vec = [vec mean(te_aFRI{elec, i})];
        end
        data_aFRI(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_aFRI(i) = ci(1);
        cihi_aFRI(i) = ci(2);
        for elec = channels_to_plot
            if (size(te_TD{elec, i}) ~= [0 0])
                vec = [vec mean(te_TD{elec, i})];
            end
        end
        if (size(vec) ~= [0 0])
            data_TD(i) = mean(vec);
            ci = bootci(1000, @mean, vec);
            cilo_TD(i) = ci(1);
            cihi_TD(i) = ci(2);
        else
            data_TD(i) = NaN;
            cilo_TD(i) = NaN;
            cihi_TD(i) = NaN;
        end
    end
    xlim(1000*[0 1.06*max(T_vec)]);
    ylim(1000*[0 0.014]);
    p0 = plot(xlim, 0.25*xlim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 3);
    p1 = plot_errorbar(1000*T_vec, 1000*data_AT, 1000*cilo_AT, 1000*cihi_AT, [1 0 0], 0, 600/40);
    p2 = plot_errorbar(1000*T_vec, 1000*data_aFRI, 1000*cilo_aFRI, 1000*cihi_aFRI, [0 0 1], 0, 600/40);
    p3 = plot_errorbar(1000*T_vec, 1000*data_TD, 1000*cilo_TD, 1000*cihi_TD, [0 1 0], 0, 600/40);
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Av. Spike Time Error (ms)', 'FontSize', 28);
    leg = legend([p0 p1 p2 p3], 'AT theoretical', 'AT', 'gAT', 'twodelta');
    set(leg, 'Location', 'NorthWest');
    set(gca, 'FontSize', 24);
    saveas(fig, [dir 't_error_mean_' name], 'epsc');
    saveas(fig, [dir 't_error_mean_' name], 'fig');

    % frac used vs T (mean + bootci)
    fig = figure;
    hold on;
    data_AT = zeros(1, length(T_vec));
    data_aFRI = zeros(1, length(T_vec));
    data_TD = zeros(1, length(T_vec));
    cihi_AT = zeros(1, length(T_vec));
    hilo_AT = zeros(1, length(T_vec));
    cihi_aFRI = zeros(1, length(T_vec));
    cilo_aFRI = zeros(1, length(T_vec));
    cihi_TD = zeros(1, length(T_vec));
    cilo_TD = zeros(1, length(T_vec));
    for i = 1 : length(T_vec)
        vec = [];
        for elec = channels_to_plot
            vec = [vec num_AT(elec, i)/den_AT(elec,i)];
        end
        data_AT(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_AT(i) = ci(1);
        cihi_AT(i) = ci(2);
        vec = [];
        for elec = channels_to_plot
            vec = [vec num_aFRI(elec, i)/den_aFRI(elec,i)];
        end
        data_aFRI(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_aFRI(i) = ci(1);
        cihi_aFRI(i) = ci(2);
        for elec = channels_to_plot
            vec = [vec num_TD(elec, i)/den_TD(elec,i)];
        end
        data_TD(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_TD(i) = ci(1);
        cihi_TD(i) = ci(2);
    end
    p1 = plot_errorbar(1000*T_vec, data_AT, cilo_AT, cihi_AT, [1 0 0], 0, 20);
    %p2 = plot_errorbar(1000*T_vec, data_aFRI, cilo_aFRI, cihi_aFRI, [0 0 1], 0, 20);
    p3 = plot_errorbar(1000*T_vec, data_TD, cilo_TD, cihi_TD, [0 1 0], 0, 20);
    leg = legend([p1 p3], 'gAT', 'twodelta', 'Location', 'Best');
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Fraction of Valid Intervals', 'FontSize', 28);
    %leg = legend([p1], 'Both methods');
    xlim(1000*[0 1.06*max(T_vec)]);
    ylim([0 1]);
    set(gca, 'FontSize', 24);
    saveas(fig, [dir 'frac_mean_' name], 'epsc');
    saveas(fig, [dir 'frac_mean_' name], 'fig');


    % frac valid for TD (1 spike, 2 spike) (mean + bootci)
    fig = figure;
    hold on;
    data_TD1 = zeros(1, length(T_vec));
    data_TD2 = zeros(1, length(T_vec));
    cihi_TD1 = zeros(1, length(T_vec));
    cilo_TD1 = zeros(1, length(T_vec));
    cihi_TD2 = zeros(1, length(T_vec));
    cilo_TD2 = zeros(1, length(T_vec));
    for i = 1:length(T_vec)
        vec1 = [];
        vec2 = [];
        for elec = channels_to_plot
            vec1 = [vec1 num_TD1(elec, i)/den_TD1(elec,i)];
            vec2 = [vec2 num_TD2(elec, i)/den_TD2(elec,i)];
        end

        data_TD1(i) = mean(vec1);
        ci = bootci(1000, @mean, vec1);
        cilo_TD1(i) = ci(1);
        cihi_TD1(i) = ci(2);

        data_TD2(i) = mean(vec2);
        if (~any(isnan(vec2)))
            ci = bootci(1000, @mean, vec2);
            cilo_TD2(i) = ci(1);
            cihi_TD2(i) = ci(2);
        else
            cilo_TD2(i) = NaN;
            cihi_TD2(i) = NaN;
        end
    end
    p1 = plot_errorbar(1000*T_vec, data_TD1, cilo_TD1, cihi_TD1, [1 0 0], 0, 20);
    p2 = plot_errorbar(1000*T_vec, data_TD2, cilo_TD2, cihi_TD2, [0 1 0], 0, 20);
    leg = legend([p1 p2], '1 spike', '2 spike', 'Location', 'Best');
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Fraction of Valid Intervals', 'FontSize', 28);
    xlim(1000*[0 1.06*max(T_vec)]);
    ylim([0 1]);
    set(gca, 'FontSize', 24);
    saveas(fig, [dir 'frac_TD_mean_' name], 'epsc');
    saveas(fig, [dir 'frac_TD_mean_' name], 'fig');




    % difference in frac used vs T (mean + bootci)
    fig = figure;
    hold on;
    data_diff = zeros(1, length(T_vec));
    cihi_diff = zeros(1, length(T_vec));
    hilo_diff = zeros(1, length(T_vec));
    for i = 1 : length(T_vec)
        vec = [];
        for elec = channels_to_plot
            vec = [vec (num_TD(elec, i)/den_TD(elec,i)-num_AT(elec,i)/den_AT(elec,i))];
        end
        data_diff(i) = mean(vec);
        ci = bootci(1000, @mean, vec);
        cilo_diff(i) = ci(1);
        cihi_diff(i) = ci(2);
    end
    p1 = plot_errorbar(1000*T_vec, data_diff, cilo_diff, cihi_diff, [1 0 0], 0, 20);
    leg = legend([p1], 'Fraction Valid of TD - Fraction Valid of AT/gAT', 'Location', 'Best');
    xlabel('Sampling Period T (ms)', 'FontSize', 28);
    ylabel('Difference in Fraction of Valid Intervals', 'FontSize', 28);
    %leg = legend([p1], 'Both methods');
    %set(leg, 'Location', 'NorthWest', 'FontSize', 28);
    xlim(1000*[0 1.06*max(T_vec)]);
    %ylim([0 1]);
    set(gca, 'FontSize', 24);
    saveas(fig, [dir 'diff_frac_mean_' name], 'epsc');
    saveas(fig, [dir 'diff_frac_mean_' name], 'fig');

end
