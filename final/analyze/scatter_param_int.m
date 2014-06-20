% Makes 3 scatter plots
% - False pos vs false neg for channels' optimal parameters
% - False pos vs false neg for channels' optimal parameters but reconstructed properly with training threshold
% - False pos vs false neg for channels with universal parameters
% And in a separate figure
% - Best params for each channel and averages

close all;
clear all;

% Parameters
method_name = 'int';
data_sets = [1, 2];

filename = [root_dir() 'analyze/opt_params.mat'];
load(filename);
obj = eval(method_name);

for data_set = data_sets

    elecs = channels_to_use(data_set);

    fp = zeros(length(elecs), 1);
    fn = zeros(length(elecs), 1);
    bestT = zeros(length(elecs), 1);
    bestL = zeros(length(elecs), 1);

    for i = 1 : length(elecs)

        elec = elecs(i);
        stats = obj{data_set, elec};
        fp(i) = stats.fp;
        fn(i) = stats.fn;
        bestT(i) = stats.bestT;
        bestL(i) = stats.bestL;

    end

    % Results with threshold estimated

    fp_est = [];
    fn_est = [];

    for i = 1 : length(elecs)
        args = struct();

        % Drew's filter
        %[A, B] = butter_hp(400, 10);
        [A, B] = butter_bp(300, 6000, 2);
        args.A_filter = A;
        args.B_filter = B;
        hp_handle = @filter_generic;  % Use causal version
        diode_handle = @ideal_diode;

        %%{
        % succ int
        T = bestT(i);
        args.L = bestL(i);
        args.K = 1;
        args.delta_t = get_dt();
        args.numiter = 3;
        args.t_resolution = T / 100;
        window = 0.005;
        %%}

        [fp_vec, fn_vec] = optthres_reconstruct(data_set, [elecs(i)], method_name, args, hp_handle, diode_handle, T, window);

        fp_est = [fp_est fp_vec];
        fn_est = [fn_est fn_vec];

    end

    % Results with threshold estimated and average universal parameters

    args = struct();

    % Drew's filter
    %[A, B] = butter_hp(400, 10);
    [A, B] = butter_bp(300, 6000, 2);
    args.A_filter = A;
    args.B_filter = B;
    hp_handle = @filter_generic;  % Use causal version
    diode_handle = @ideal_diode;

    %%{
    % succ int
    T = mean(bestT);
    args.L = round(mean(bestL));
    args.K = 1;
    args.delta_t = get_dt();
    args.numiter = 3;
    args.t_resolution = T / 100;
    window = 0.005;
    %%}

    [fp_uni, fn_uni] = optthres_reconstruct(data_set, elecs, method_name, args, hp_handle, diode_handle, T, window);

    fig = figure;
    subplot(2, 1, 1);
    hold on;
    for i = 1 : length(fn)
        plot([fn(i), fn_est(i), fn_uni(i)], [fp(i), fp_est(i), fp_uni(i)], 'k:');
    end
    p1 = plot(fn, fp, 'b.');
    p2 = plot(fn_est, fp_est, 'r.');
    p3 = plot(fn_uni, fp_uni, 'g.');
    xlabel('false negatives');
    ylabel('false positives');
    title(['data set ' int2str(data_set)]);
    xlim([0, 1]);
    ylim([0, 1]);
    legend([p1, p2, p3], 'best params, best threshold', 'best params, estimated threshold', 'universal params, estimated threshold', 'Location', 'EastOutside');
    axis square;

    subplot(2, 1, 2);
    hold on;
    plot(bestT, bestL, 'b.');
    plot(mean(bestT), mean(bestL), 'r.');
    xlabel('T');
    ylabel('L');
    title(['optimal params: T = ' num2str(mean(bestT)) ', L = ' num2str(mean(bestL))]);
    legend('channel-by-channel', 'average', 'Location', 'EastOutside');
    axis square;

    saveas(fig, ['~/Desktop/scatter_' int2str(data_set) '.pdf']);

end
