% Script for generating RMSE vs sampling rate plots

clear all;
close all;

%{

for i = 20:22
    reconstruct_and_save('ann', i, [1]);
end

%}

%{

for i = 15 : 16
    reconstruct_and_save('ann', i, [1]);
end

%}

%reconstruct_and_save('ann', 1, [1]);

%{

clear all;

for i = 50 : 59
    reconstruct_and_save('gibbs', i, [1]);
end

%}

%{

% Richardson figure, Gibbs

data_set = [2, 5:17, 19];
for i = data_set
    reconstruct_and_save('gibbs', 34, i, 1);
end

%}

%est_delay_b('ann', 1, 1);

%{

paramIDs = 10 : 16;
method_names = {'gibbs'};  %{'gibbs', 'RMSE', 'ker'};

for i = 1 : length(method_names)
    method_name = method_names{i};
    est_delay_b(method_name, paramIDs, 1);
    close all;
end

%}

%{

paramIDs = 20 : 33;
method_names = {'RMSE'};

for i = 1 : length(method_names)
    method_name = method_names{i};
    est_delay_b(method_name, paramIDs, 1);
    close all;
end

paramIDs = 20 : 36;
method_names = {'gibbs'};

for i = 1 : length(method_names)
    method_name = method_names{i};
    est_delay_b(method_name, paramIDs, 1);
    close all;
end

%}

%{

% Compare 4 methods plot

% NOTE: this should match uni_params.m
sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000];%, 3000, 4000, 5000];
paramIDs = 10 : 16;

method_names = {'ker', 'RMSE', 'gibbs', 'ML', 'RMSE_cal'};

to_plot = [];  % stuff to plot
ci_to_plot = [];  % confidence intervals to plot

for i = 1 : length(method_names)
    
    method_name = method_names{i};

    RMSE = zeros(2, 0);
    ci_r = zeros(2, 0); % real
    ci = zeros(2, 0); % reconstructed
    for paramID = paramIDs
        [RMSE_real, RMSE_rec, ci_real, ci_rec] = run_decode(method_name, paramID);
        RMSE_real
        RMSE_rec
        RMSE = [RMSE [RMSE_real ; RMSE_rec]];
        RMSE
        ci_r = [ci_r ci_real'];
        ci = [ci ci_rec'];
        %close all;
    end
    
    % put ci's in terms of difference from mean
    ci_r = (ci_r(1, :) - ci_r(2, :)) / 2;
    ci = (ci(1, :) - ci(2, :)) / 2;
    
    if i == 1
        to_plot = RMSE;
        ci_to_plot = [ci_r ; ci];
    else
        to_plot = [to_plot ; RMSE(2, :)];
        ci_to_plot = [ci_to_plot ; ci];
    end
    %plot(sampling_rates, ci_r, 'c');  % real ci
    %plot(sampling_rates, ci, 'y');  % rec ci
    
end

% annihilating filter
[RMSE_real, RMSE_rec, ci_real, ci_rec] = run_decode('ann', 1);
ann_RMSE = RMSE_rec

close all;
figure;
hold on;
errorbar(sampling_rates, to_plot(1, :), ci_to_plot(1, :), 'b');
errorbar(sampling_rates, to_plot(2, :), ci_to_plot(2, :), 'r');
errorbar(sampling_rates, to_plot(3, :), ci_to_plot(3, :), 'g');
errorbar(sampling_rates, to_plot(4, :), ci_to_plot(4, :), 'c');
errorbar(sampling_rates, to_plot(5, :), ci_to_plot(5, :), 'm');
errorbar(sampling_rates, to_plot(6, :), ci_to_plot(6, :), 'k');
xlabel('sampling rate');
ylabel('RMSE');
legend('real', 'ker', 'RMSE', 'gibbs', 'ML', 'RMSE_cal');
y_lim = ylim;
ylim([0, y_lim(2)]);

%}

%{

% RMSE rectification comparison plot

% NOTE: this should match uni_params.m
sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000];%, 3000, 4000, 5000];
paramIDs = 10 : 16;

method_names = {'RMSE', 'RMSE'};

to_plot = [];  % stuff to plot
ci_to_plot = [];  % confidence intervals to plot

for i = 1 : length(method_names)
    
    method_name = method_names{i};

    RMSE = zeros(2, 0);
    ci_r = zeros(2, 0); % real
    ci = zeros(2, 0); % reconstructed
    for paramID = paramIDs
        p_arg = paramID;
        if i == 1  % do unrectified first
            p_arg = p_arg + 10;
        end
        [RMSE_real, RMSE_rec, ci_real, ci_rec] = run_decode(method_name, p_arg);
        RMSE_real
        RMSE_rec
        RMSE = [RMSE [RMSE_real ; RMSE_rec]];
        RMSE
        ci_r = [ci_r ci_real'];
        ci = [ci ci_rec'];
        %close all;
    end
    
    % put ci's in terms of difference from mean
    ci_r = (ci_r(1, :) - ci_r(2, :)) / 2;
    ci = (ci(1, :) - ci(2, :)) / 2;
    
    if i == 1
        to_plot = RMSE;
        ci_to_plot = [ci_r ; ci];
    else
        to_plot = [to_plot ; RMSE(2, :)];
        ci_to_plot = [ci_to_plot ; ci];
    end
    %plot(sampling_rates, ci_r, 'c');  % real ci
    %plot(sampling_rates, ci, 'y');  % rec ci
    
end

figure;
hold on;
errorbar(sampling_rates, to_plot(1, :), ci_to_plot(1, :), 'b');
errorbar(sampling_rates, to_plot(2, :), ci_to_plot(2, :), 'r');
errorbar(sampling_rates, to_plot(3, :), ci_to_plot(3, :), 'g');
%errorbar(sampling_rates, to_plot(4, :), ci_to_plot(4, :), 'c');
xlabel('sampling rate');
ylabel('RMSE');
legend('real', 'RMSE', 'RMSE with half-wave rectification');
y_lim = ylim;
ylim([0, y_lim(2)]);

%}

%%{

% Gibbs plot for conference

% NOTE: this should match uni_params.m
sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000];%, 3000, 4000, 5000];
paramIDs = 0 : 6;

method_names = {'ann', 'gibbs'};%, 'RMSE'};
paramOffset = [10, 30, 10];

to_plot = [];  % stuff to plot
ci_to_plot = [];  % confidence intervals to plot

for i = 1 : length(method_names)
    
    method_name = method_names{i};

    RMSE = zeros(2, 0);
    ci_r = zeros(2, 0); % real
    ci = zeros(2, 0); % reconstructed
    ci_s = zeros(2, 0); % still
    data_still = [];
    for paramID = paramIDs + paramOffset(i)
        
        [RMSE_real, RMSE_rec, ci_real, ci_rec, RMSE_still, ci_still] = run_decode(method_name, paramID);
        
        RMSE_real
        RMSE_rec
        RMSE = [RMSE [RMSE_real ; RMSE_rec]];
        RMSE
        ci_r = [ci_r ci_real'];
        ci = [ci ci_rec'];
        ci_s = [ci_s ci_still'];
        %close all;
        data_still = [data_still RMSE_still];
    end
    
    if i == 1
        to_plot = RMSE;
        ci_to_plot = [ci_r ; ci];
    else
        to_plot = [to_plot ; RMSE(2, :)];
        ci_to_plot = [ci_to_plot ; ci];
    end
    %plot(sampling_rates, ci_r, 'c');  % real ci
    %plot(sampling_rates, ci, 'y');  % rec ci
    
end

save('C:\Users\alex\Desktop\figures\mat\rmse');

%%}

%%{

load('C:\Users\alex\Desktop\figures\mat\rmse');

close all;
f = figure;
hold on;
p1 = plot_errorbar(sampling_rates, to_plot(1, :), ci_to_plot(2, :), ci_to_plot(1, :), [0 0 0]);
p2 = plot_errorbar(sampling_rates, to_plot(2, :), ci_to_plot(4, :), ci_to_plot(3, :), [0 0 1]);
p3 = plot_errorbar(sampling_rates, to_plot(3, :), ci_to_plot(6, :), ci_to_plot(5, :), [0 1 0]);
%p4 = plot_errorbar(sampling_rates, to_plot(4, :), ci_to_plot(8, :), ci_to_plot(7, :), [1 0 1]);
p5 = plot_errorbar(sampling_rates, data_still, ci_s(1, :), ci_s(2, :), [1 0 0]);
%errorbar(sampling_rates, to_plot(3, :), ci_to_plot(3, :), 'g');
%errorbar(sampling_rates, to_plot(4, :), ci_to_plot(4, :), 'c');
xlabel('sampling rate (Hz)', 'FontSize', 18);
ylabel('RMSE in Position (cm)', 'FontSize', 18);
leg = legend([p1 p2 p3 p5], 'real spikes', 'annihilating filter', 'gibbs', 'no motion');
set(leg, 'FontSize', 18, 'Location', 'SouthEast');
y_lim = ylim;
ylim([0, y_lim(2)]);

saveas(f, ['C:\Users\alex\Desktop\rmse_sampl'], 'epsc');

%%}

%{

for i = 20:22
    [RMSE_real, RMSE_rec, ci_real, ci_rec, RMSE_still, ci_still] = run_decode('ann', i);
    RMSE_rec
end

%}