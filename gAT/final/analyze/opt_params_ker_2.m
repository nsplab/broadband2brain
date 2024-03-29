% Optimizes parameters and saves best params in .mat file
% Uses fixed T, plots performance against T/L and kernel width

close all;
clear all;

% Format for .mat file
% Cell array of name method_name
% method_name{data_set, electrode} = struct with data
filename = [root_dir() 'analyze/opt_params.mat'];
load(filename);

tic

% Parameters
data_sets = [1, 2];
electrodes = [16]; %channels_to_use(data_set);

%{
args = struct();

% Drew's filter
[A, B] = butter_hp(1000, 2);
%[A, B] = butter_bp(300, 6000, 2);
args.A_filter = A;
args.B_filter = B;
hp_handle = @filter_generic;  % Use causal version
diode_handle = @ideal_diode;
%}

%{
% succ int
T = 0.02;
args.L = 3;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'int';
window = 0.005;
%}

%{
% sinc
T = 0.02;
args.B = 300;
args.L = 10;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'ker';
args.h = @normalized_sinc;
window = 0.1;
%}

%{
% RC
T = 0.02;
args.RC = 0.003;
args.L = 10;
args.K = 1;
args.delta_t = get_dt();
args.numiter = 3;
args.t_resolution = T / 100;
method_name = 'ker';
args.h = @RC_lowpass;
window = 0.1;
%}

method_name = 'ker';
paramID = 2;
[args, T, window, hp_handle, diode_handle, t0] = uni_params(method_name, paramID);

%sampling_rates = [200, 400, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000];
%sampling_rates = 200 : 200 : 5000;
sampling_rates = 50 : 50 : 1500;
%sampling_rates = [200, 400, 600];
Ls = 1 ./ sampling_rates;  % T_s
T_s_vec = Ls
%Ls = linspace(T / 120, T / 5, 8);  % (should be called 'T_s')  TODO: choose
%Ts = 0.01 : 0.0005 : 0.04;
%RCs = 0.001 : 0.001 : 0.01;
Ts = linspace(1/1000, 25/1000, 25);  % width

for data_set = data_sets

    % Do all electrodes
    electrodes = channels_to_use(data_set);

    for elec = electrodes

        % Optimize RC

        %{
        err_rate = zeros(length(RCs), 1);

        for i = 1 : length(RCs)

            disp([int2str(i) ' of ' int2str(length(RCs))]);

            args.RC = RCs(i);

            [true_pos_rate, false_pos_rate, false_neg_rate] = compare_to_real(data_set, elec, hp_handle, diode_handle, method_name, args, T, window);

            error_rate = false_pos_rate + false_neg_rate;

            err_rate(i) = error_rate;

        end

        figure;
        plot(RCs, err_rate, 'b');
        xlabel('RC');
        ylabel('error rate');
        title(['method: ' method_name ', channel: ' int2str(elec)]);
        %}

        % Heat plot of T and L
        count = 0;
        err_rate = zeros(length(Ls), length(Ts));
        fp_rate = zeros(length(Ls), length(Ts));        
        fn_rate = zeros(length(Ls), length(Ts));
        for i = 1 : length(Ls)
            for j = 1 : length(Ts)
                count = count + 1;
                disp([int2str(count) ' of ' int2str(length(Ls)*length(Ts))]);
                args.L = round(T / Ls(i));
                %T = Ls(i) * args.L;
                %T = Ts(j);
                args.s1 = -5 / Ts(j);
                args.s2 = -4 / Ts(j);
                [true_pos_rate, false_pos_rate, false_neg_rate] = compare_to_real(data_set, elec, hp_handle, diode_handle, method_name, args, T, window);
                fp_rate(i, j) = false_pos_rate;
                fn_rate(i, j) = false_neg_rate;
            end
        end

        err_rate = fp_rate + fn_rate;

        % Find min
        [val, ind_L] = min(err_rate);
        [min_err, ind_T] = min(val);
        ind_L = ind_L(ind_T);
        bestT = Ts(ind_T);
        bestL = sampling_rates(ind_L); %Ls(ind_L);
        fp = fp_rate(ind_L, ind_T);
        fn = fn_rate(ind_L, ind_T);
        
        title_str = {['method: ' method_name ', data set: ' int2str(data_set) ', channel: ' int2str(elec)], ...
            ['min err. rate = ' num2str(min_err) ', false pos. rate = ' num2str(fp) ', false neg. rate = ' num2str(fn)], ...
            ['best width = ' num2str(bestT) ', best T_s = ' num2str(bestL)]};

        % Heat map
        fig = figure;
        subplot(2, 2, 1);
        hold on;
        image(Ts, sampling_rates, err_rate * length(get(gcf, 'Colormap')));
        axis tight;
        plot([bestT, bestT], ylim, 'm--');
        plot(xlim, [bestL, bestL], 'g--');
        plot([bestT], [bestL], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        set(gca,'ydir','normal')        
        xlabel('width');
        ylabel('sampl. rate');
        zlabel('error rate');
        title(title_str);

        % 3D surface plot
        subplot(2, 2, 2);
        surf(Ts, sampling_rates, err_rate);
        caxis([0, 1]);
        xlabel('width');
        ylabel('sampl. rate');
        zlabel('error rate');
        %title(title_str);
        axis tight;
        zlim([0, 1]);

        % Projections
        subplot(2, 2, 3);
        hold on;
        plot(Ts, err_rate(ind_L, :), 'g');
        xlabel('width');
        ylabel('error rate');
        plot([bestT], [min_err], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

        subplot(2, 2, 4);
        hold on;
        plot(sampling_rates, err_rate(:, ind_T), 'm');
        xlabel('sampl. rate');
        ylabel('error rate');
        plot([bestL], [min_err], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

        % Save
        obj = [method_name '{' int2str(data_set) ', ' int2str(elec) '}'];
        eval([obj ' = struct();']);
        eval([obj '.Ts = Ts;']);
        eval([obj '.Ls = Ls;']);
        eval([obj '.err_rate = err_rate;']);
        eval([obj '.fp_rate = fp_rate;']);
        eval([obj '.fn_rate = fn_rate;']);
        eval([obj '.min_err = min_err;']);
        eval([obj '.fp = fp;']);
        eval([obj '.fn = fn;']);
        eval([obj '.bestT = bestT;']);
        eval([obj '.bestL = bestL;']);
        eval([obj '.ind_T = ind_T;']);
        eval([obj '.ind_L = ind_L;']);

        % Save figure
        mkdir([root_dir() 'opt_params_ker_2']);
        saveas(fig, [root_dir() 'opt_params_ker_2/' method_name '_' int2str(data_set) '_' int2str(elec) '.pdf']);

        % Save results
        save(filename, method_name, '-append');
        
        close;

    end
end

toc
