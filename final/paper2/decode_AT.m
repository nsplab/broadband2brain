% RMES vs sampling rate for AT, aFRI, real, no motion

% NOTE: this should match uni_params.m
T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];
paramIDs = 10:20;

% exclude first 3
T_vec = T_vec(4:end);
paramIDs = paramIDs(4:end);

sampling_rates = 1./T_vec;

method_names = {'AT', 'analog'};
paramOffset = [0 0];

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

save('C:\Users\alex\Desktop\figures\mat\rmse_AT');

%%}

%%{

load('C:\Users\alex\Desktop\figures\mat\rmse_AT');

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
leg = legend([p1 p2 p3 p5], 'real spikes', 'AT', 'aFRI', 'no motion');
set(leg, 'FontSize', 18, 'Location', 'SouthEast');
y_lim = ylim;
ylim([0, y_lim(2)]);

%saveas(f, ['C:\Users\alex\Desktop\rmse_sampl_AT'], 'epsc');