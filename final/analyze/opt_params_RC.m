% Optimizes parameters and saves best params in .mat file

close all;
clear all;

% Format for .mat file
% Cell array of name method_name
% method_name{data_set, electrode} = struct with data
filename = [root_dir() 'analyze/opt_params.mat'];
load(filename);

tic

% Parameters
sampling_rates = [500, 1000, 1500];
data_sets = [1, 2];
electrodes = [9]; %channels_to_use(data_set);

args = struct();

% Drew's filter
%[A, B] = butter_hp(400, 10);
[A, B] = butter_bp(300, 6000, 2);
args.A_filter = A;
args.B_filter = B;
hp_handle = @filter_generic;  % Use causal version
diode_handle = @ideal_diode;

%{
% succ int
T = 0.02;
args.L = 8;
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

%%{
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
%%}

%Ls = 2 : 1 : 15;
%Ts = 0.01 : 0.002 : 0.04;
RCs = 0.001 : 0.001 : 0.01;
min_T = 0.01;
max_T = 0.04;

for sr_ind = 1 : length(sampling_rates)

    sampling_rate = sampling_rates(sr_ind);

    Ls = ceil(sampling_rate * min_T) : round(sampling_rate / 500) : floor(sampling_rate * max_T);

    Ts = Ls / sampling_rate;

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

            % Heat plot of T and RC
            count = 0;
            err_rate = zeros(length(Ls), length(RCs));
            fp_rate = zeros(length(Ls), length(RCs));        
            fn_rate = zeros(length(Ls), length(RCs));
            for i = 1 : length(Ls)
                for j = 1 : length(RCs)
                    count = count + 1;
                    disp([int2str(count) ' of ' int2str(length(Ls)*length(RCs))]);
                    args.L = Ls(i);
                    T = Ts(i);
                    args.RC = RCs(j);
                    [true_pos_rate, false_pos_rate, false_neg_rate] = compare_to_real(data_set, elec, hp_handle, diode_handle, method_name, args, T, window);
                    fp_rate(i, j) = false_pos_rate;
                    fn_rate(i, j) = false_neg_rate;
                end
            end

            err_rate = fp_rate + fn_rate;

            % Find min
            [val, ind_T] = min(err_rate);
            [min_err, ind_RC] = min(val);
            ind_T = ind_T(ind_RC);
            bestT = Ts(ind_T);
            bestL = Ls(ind_T);
            bestRC = RCs(ind_RC);
            fp = fp_rate(ind_T, ind_RC);
            fn = fn_rate(ind_T, ind_RC);

            title_str = {['method: ' method_name ', data set: ' int2str(data_set) ', channel: ' int2str(elec)], ...
                ['min err. rate = ' num2str(min_err) ', false pos. rate = ' num2str(fp) ', false neg. rate = ' num2str(fn)], ...
                ['best T = ' num2str(bestT) ', best RC = ' num2str(bestRC) ', sampling rate = ' num2str(bestL/bestT) ' Hz']};

            % Heat map
            fig = figure;
            subplot(2, 2, 1);
            hold on;
            image(RCs, Ts, err_rate * length(get(gcf, 'Colormap')));
            axis tight;
            plot([bestRC, bestRC], ylim, 'm--');
            plot(xlim, [bestT, bestT], 'g--');
            plot([bestRC], [bestT], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
            set(gca,'ydir','normal')        
            xlabel('RC');
            ylabel('T');
            zlabel('error rate');
            title(title_str);

            % 3D surface plot
            subplot(2, 2, 2);
            surf(RCs, Ts, err_rate);
            caxis([0, 1]);
            xlabel('RC');
            ylabel('T');
            zlabel('error rate');
            %title(title_str);
            axis tight;
            zlim([0, 1]);

            % Projections
            subplot(2, 2, 3);
            hold on;
            plot(RCs, err_rate(ind_T, :), 'g');
            xlabel('RC');
            ylabel('error rate');
            plot([bestRC], [min_err], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

            subplot(2, 2, 4);
            hold on;
            plot(Ts, err_rate(:, ind_RC), 'm');
            xlabel('T');
            ylabel('error rate');
            plot([bestT], [min_err], 'Marker', 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

            % Save
            obj = [method_name '{' int2str(data_set) ', ' int2str(elec) ', ' int2str(sr_ind) '}'];
            eval([obj ' = struct();']);
            eval([obj '.Ts = Ts;']);
            eval([obj '.Ls = Ls;']);
            eval([obj '.RCs = RCs;']);
            eval([obj '.err_rate = err_rate;']);
            eval([obj '.fp_rate = fp_rate;']);
            eval([obj '.fn_rate = fn_rate;']);
            eval([obj '.min_err = min_err;']);
            eval([obj '.fp = fp;']);
            eval([obj '.fn = fn;']);
            eval([obj '.bestT = bestT;']);
            eval([obj '.bestL = bestL;']);
            eval([obj '.bestRC = bestRC;']);
            eval([obj '.ind_T = ind_T;']);
            eval([obj '.ind_RC = ind_RC;']);

            % Save figure
            saveas(fig, ['~/Desktop/opt_ker/opt_' method_name '_' int2str(data_set) '_' int2str(elec) '_' num2str(sampling_rate) '.pdf']);

            % Save results
            save(filename, method_name, '-append');

        end
    end
end

toc
