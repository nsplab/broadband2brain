% Compares AT and A-FRI on a spike-by-spike basis
% Plots # failures (false positives + false negatives) and av_t_error vs
% sampling period
function [] = spk_by_spk(channel)

close all;

% Parameters
data_set = 1;
electrode = 16;
if nargin >= 1
    electrode = channel
end
duration = 60;
delta = 0.0011;  % refractory period

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
av_spk_rate = length(real_spk)/duration

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % Sampling periods
num_thres = 20;
thres_vec = linspace(0.9*mean(dat)+0.1*max(dat), 0.3*mean(dat)+0.7*max(dat), num_thres);  % thresholds to try

% Loop through methods
methods = {'Analog Thresholding', 'aFRI'};
colors = {[1 0 0], [0 0 1]};

% 2 overall figures
fail_fig = figure;
hold on;
xlabel('sampling period (T)');
ylabel('fail rate');

t_fig = figure;
hold on;
xlabel('sampling period (T)');
ylabel('av. t error');
yx4 = plot(T_vec, T_vec/4, 'k--');

fpfn_figs = zeros(1, length(T_vec));
for i = 1 : length(fpfn_figs)
    fpfn_figs(i) = figure;
    hold on;
    xlabel('false negatives per true spike');
    ylabel('false positives per true spike');
    ylim([0 1]);
end

roc_figs = zeros(1, length(T_vec));
for i = 1 : length(roc_figs)
    roc_figs(i) = figure;
    hold on;
    xlabel('false negatives per true spike');
    ylabel('false positives per true spike');
    ylim([0 1]);
end

h_vec = [yx4];

for m = 1 : length(methods)

    fail_vec = zeros(1, length(T_vec));
    av_delta_t_vec = zeros(1, length(T_vec));

    % Loop through T
    for i = 1 : length(T_vec)
        T = T_vec(i);
        t_block = start_time : T : end_time;
        t_block = t_block(1 : end-1);  % t_block contains start times
        
        fp_vec = zeros(1, length(thres_vec));
        fn_vec = zeros(1, length(thres_vec));
        delta_t_vec = zeros(1, length(thres_vec));
        
        fp_roc = zeros(1, length(thres_vec));
        fn_roc = zeros(1, length(thres_vec));
        
        % Loop through thres
        for j = 1 : length(thres_vec)
            thres = thres_vec(j);

            delta_t_sum = 0;
            delta_t_count = 0;
            
            real_ind = 1;  % index in real_spk
            
            % TESTING
            plotflag = 0;
            %{
            target_thres = 75.5388;
            plotflag = (T == 0.02 && abs(thres - target_thres) < 0.001 && m == 2);
            if plotflag
                testfig = figure;
                hold on;
                plot(time, dat, 'c');
                plot(t_block, ones(size(t_block)), 'b.');
                plot(real_spk, -5*ones(size(real_spk)), 'k.');
            end
            %}
            
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
                
                % get reconstructed spikes in interval
                % analog thresholding
                t_AT = T/2;
                num_rec_spikes_AT = max(seg) > thres;
                % aFRI
                seg_comp = seg > thres;
                y = succ_int(seg_comp, get_dt(), 2);
                [t_F, c] = reconstruct_analog(y, T, [], []);
                num_rec_spikes_F = (c ~= 0);
                if m == 1
                    t = t_AT;
                    num_rec_spikes = num_rec_spikes_AT;
                elseif m == 2
                    t = t_F;
                    num_rec_spikes = num_rec_spikes_F;
                end
                
                if num_rec_spikes < num_real_spikes
                    fn_vec(j) = fn_vec(j) + num_real_spikes - num_rec_spikes;
                elseif num_rec_spikes > num_real_spikes
                    fp_vec(j) = fp_vec(j) + num_rec_spikes - num_real_spikes;
                elseif num_rec_spikes == 1
                    delta_t_count = delta_t_count + 1;
                    delta_t_sum = delta_t_sum + abs(t - t_real);
                end
                
                if num_rec_spikes == 1 && (length(t_rec) < 1 || t_block(k)+t > t_rec(end)+delta)
                    t_rec = [t_rec t_block(k)+t];
                end
                
                % TESTING
                if plotflag && num_rec_spikes == 1
                    figure(testfig);
                    plot(t + t_block(k), -10, 'm.');
                end
                
            end
            
            delta_t_vec(j) = delta_t_sum / delta_t_count;
            
            % for ROC
            tolerance = 0.0025;
            [~, ~, ~, ~, ~, fp_rate, fn_rate] = compare_spikes(t_rec, real_spk, tolerance);
            fp_roc(j) = fp_rate;
            fn_roc(j) = fn_rate;
            
        end
        
        % For this T value, find best threshold
        [val ind] = min(fp_vec + fn_vec);
        fail_vec(i) = val;
        av_delta_t_vec(i) = delta_t_vec(ind);

        % Plot fp/fn curve for this T value
        figure(fpfn_figs(i));
        plot_errorbar(fn_vec/length(real_spk), fp_vec/length(real_spk), fp_vec/length(real_spk), fp_vec/length(real_spk), colors{m}, 1);
        plot_dots(fn_vec(ind)/length(real_spk), fp_vec(ind)/length(real_spk), [0 1 0]);
        title(['T = ' num2str(T) ', best thres = ' num2str(thres_vec(ind)) ', error rate = ' num2str(val/length(real_spk))]);
        
        % Plot ROC curve using raw compare_spikes
        [val2 ind2] = min(fp_roc + fn_roc);
        figure(roc_figs(i));
        plot_errorbar(fn_roc, fp_roc, fp_roc, fp_roc, colors{m}, 1);
        plot_dots(fn_roc(ind2), fp_roc(ind2), [0 1 0]);
        title(['ROC: T = ' num2str(T)]);
        %title(['ROC: T = ' num2str(T) ', best thres = ' num2str(thres_vec(ind)) ', error rate = ' num2str(val/length(real_spk))]);
        p1 = plot_errorbar(-10, -10, -10, -10, colors{1}, 1);
        p2 = plot_errorbar(-10, -10, -10, -10, colors{2}, 1);
        legend([p1 p2], 'AT', 'aFRI');
        
    end
    
    % Plot
    figure(fail_fig);
    p = plot_errorbar(T_vec, fail_vec/length(real_spk), fail_vec/length(real_spk), fail_vec/length(real_spk), colors{m}, 1);
    legend(p, 'both methods');
    figure(t_fig);
    p1 = plot_errorbar(T_vec, av_delta_t_vec, av_delta_t_vec, av_delta_t_vec, colors{m}, 1);
    h_vec = [h_vec p1];
    legend(h_vec, 'y=x/4', 'Analog Thresholding', 'aFRI');
    
end