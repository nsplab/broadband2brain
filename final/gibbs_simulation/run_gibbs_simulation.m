% Fair comparison in simulated environment between Gibbs sampling and iterative ML
% Noise added to samples not to spikes

function [] = run_gibbs_simulation(mode, compute)

close all;

% Parameters
T = 20;  % [0, T] is time window containing spikes
K = 2;  % number of spikes
N = 10;  % number of samples
sigmah = 2;
sigmae = 1;  %2;
numiter = 50;% 50 or 3
mean_iter = 25; % 25 or 2
window = 10;

dir = 'C:\Users\alex\Desktop\figures\';

% CHANGE MODE HERE
if nargin < 1
    mode = 0;
    compute = 0;  % 1 to compute, 0 to load from .mat files
end

if mode == 0
    % Single run with default params, debug mode (plots)
    [error_gibbs error_ml] = gibbs_simulation(T, K, N, sigmah, sigmae, numiter, mean_iter, window, 1);
elseif mode == 1
    % Average over a bunch of runs with default params
    e_g_av = 0;
    e_ga_av = 0;
    e_ml_av = 0;
    like_av = 0;
    err_av = 0;
    l4_av = 0;
    l5_av = 0;
    numtrials = 100;
    for i = 1 : numtrials
        [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2] = gibbs_simulation(T, K, N, sigmah, sigmae, numiter, mean_iter, window, 0);
        e_g_av = e_g_av + error_g;
        e_ga_av = e_ga_av + error_ga;
        e_ml_av = e_ml_av + error_ml;
        like_av = like_av + flag1;
        err_av = err_av + flag2;
        %l4_av = l4_av + LIKE_gibbv4;
        %l5_av = l5_av + LIKE_gibbv5;
    end
    e_g_av = e_g_av / numtrials
    e_ga_av = e_ga_av / numtrials
    e_ml_av = e_ml_av / numtrials
    like_av = like_av / numtrials
    err_av = err_av / numtrials
    %like_gibbv4_av = l4_av / numtrials
    %like_gibbv5_av = l5_av / numtrials
elseif mode == 2
    % Error vs numiter
    if compute
        itervec = [1:9 10:10:100];
        numtrials = 100;
        e_g_vec = zeros(size(itervec));
        e_ga_vec = zeros(size(itervec));
        e_ml_vec = zeros(size(itervec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(itervec)
            disp([int2str(i) ' of ' int2str(length(itervec))]);
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                mean_iter = ceil(itervec(i) / 2);  % CHANGE MEAN_ITER HERE
                [error_gibbs error_ga error_ml] = gibbs_simulation(T, K, N, sigmah, sigmae, itervec(i), mean_iter, window, 0);
                e_g(j) = error_gibbs;
                e_ga(j) = error_ga;
                e_ml(j) = error_ml;
            end
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);%[prctile(e_g, 5) ; prctile(e_g, 95)];
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    % Prune
    ind = [1 2 3 5 7 10:length(itervec)];
    itervec = itervec(ind);
    e_g_vec = e_g_vec(ind);
    e_ga_vec = e_ga_vec(ind);
    e_ml_vec = e_ml_vec(ind);
    e_g_ci = e_g_ci(:,ind);
    e_ga_ci = e_ga_ci(:,ind);
    e_ml_ci = e_ml_ci(:,ind);
    % Full
    x1 = -2;
    x2 = 54;
    y1 = 0.002;
    y2 = 0.098;
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(itervec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs',1);
    p2 = plot_errorbar(itervec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga',1);
    p3 = plot_errorbar(itervec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml',1);
    plot_line([x1 x2 x2 x1 x1], [y2 y2 y1 y1 y2], 'err', 1, 2);
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('No. Iterations ($I$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    xlim([-3 105]);
    ylim([0 0.55]);
    saveas(f, [dir 'iter_full'], 'epsc');
    % Zoom
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(itervec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');%, 1);
    p2 = plot_errorbar(itervec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');%, 1);
    p3 = plot_errorbar(itervec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');%, 1);
    plot_line([x1 x2 x2 x1 x1], [y2 y2 y1 y1 y2], 'err', 1, 2);
    %leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    %set(leg, 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('No. Iterations ($I$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    xlim([-3 55]);
    ylim([0 0.1]);
    saveas(f, [dir 'iter_zoom'], 'epsc');
elseif mode == 3
    % Error vs sigmae
    if compute
        sigvec = [0.001 0.5:0.5:8];
        numtrials = 100;
        e_g_vec = zeros(size(sigvec));
        e_ga_vec = zeros(size(sigvec));
        e_ml_vec = zeros(size(sigvec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        gibbs_all = zeros(0, numtrials);
        for i = 1 : length(sigvec)
            i
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_gibbs error_ga error_ml] = gibbs_simulation(T, K, N, sigmah, sigvec(i), numiter, mean_iter, window, 0);
                e_g(j) = error_gibbs;
                e_ga(j) = error_ga;
                e_ml(j) = error_ml;
            end
            gibbs_all = [gibbs_all ; e_g];
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(sigvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs',1);
    p2 = plot_errorbar(sigvec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga',1);
    p3 = plot_errorbar(sigvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml',1);
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    %title('Error vs Noise Level', 'FontSize', 14);
    xlim([-0.5, 8.5]);
    ylim([-0.05, 0.8]);
    saveas(f, [dir 'noise'], 'epsc');
elseif mode == 4
    % Error vs N (# samples)
    if compute
        Nvec = [10:5:20 30:10:100];
        numtrials = 100;
        e_g_vec = zeros(size(Nvec));
        e_ga_vec = zeros(size(Nvec));
        e_ml_vec = zeros(size(Nvec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(Nvec)
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_gibbs error_ga error_ml] = gibbs_simulation(T, K, Nvec(i), sigmah, sigmae, numiter, mean_iter, window, 0);
                e_g(j) = error_gibbs;
                e_ga(j) = error_ga;
                e_ml(j) = error_ml;
            end
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(Nvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(Nvec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(Nvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('No. Samples ($N$)', 'Interpreter', 'latex', 'FontSize', 24);% 24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    ylim([0 0.35]);
    xlim([0 105]);
    %title('Error vs Number of Samples', 'FontSize', 14);
    saveas(f, [dir 'samples'], 'epsc');
elseif mode == 5
    % Evaluate estimation of sigmae
    % Plot sigmae_hat (by both methods) vs sigmae
    % Error vs sigmae
    sig_mode = 1;  % 0: std, 1: E[theta | y]
    if compute
        sigvec = [0.001 0.5:0.5:8];
        numtrials = 1000;
        e_g_vec = zeros(size(sigvec));
        e_ml_vec = zeros(size(sigvec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(sigvec)
            i
            disp([num2str(i) ' of ' num2str(length(sigvec))]);
            e_g = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_g error_ga error_ml sigma_g sigma_ml] = gibbs_simulation(T, K, N, sigmah, sigvec(i), numiter, mean_iter, window, 0, sig_mode);
                e_g(j) = sigma_g;
                e_ml(j) = sigma_ml;
            end
            
            % cap sigma
            sig_max = 25;
            e_g = min(e_g, sig_max);
            e_ml = min(e_ml, sig_max);
            
            e_g_vec(i) = mean(e_g);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end

        factor_vec = e_ml_vec ./ sigvec
        mean(factor_vec(2:end))
        
        save([dir 'mat/save' int2str(mode) int2str(sig_mode)]);
    else
        load([dir 'mat/save' int2str(mode) '0']);
        e_ml_vec_0 = e_ml_vec;
        e_ml_ci_0 = e_ml_ci;
        load([dir 'mat/save' int2str(mode) '1']);
        e_ml_vec_1 = e_ml_vec;
        e_ml_ci_1 = e_ml_ci;
    end
    
    f = figure;
    box on;
    hold on;
    p3 = plot_errorbar(sigvec, e_ml_vec_1, e_ml_ci_1(1,:), e_ml_ci_1(2,:), lighten('iterml'), 1);
    p1 = plot_errorbar(sigvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs', 1);
    p2 = plot_errorbar(sigvec, e_ml_vec_0, e_ml_ci_0(1,:), e_ml_ci_0(2,:), string_to_color('iterml')*0.9, 1);
    h = plot([0, max(sigvec)], [0, max(sigvec)], 'k:');
    set(h, 'LineWidth', 5);
    set(gca, 'FontSize', 18);
    leg = legend([p1 p2 p3], 'Gibbs', 'IterML (I)', 'IterML (II)');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('Estimated Noise Level ($\widehat{\sigma_e}$)', 'Interpreter', 'latex', 'FontSize', 24);
    %if sig_mode == 0
    %    title('$\widehat{\sigma_e} = \mathrm{stdev}({\bf y}-{\bf z})$', 'FontSize', 28, 'Interpreter', 'latex');
    %else
    %    title('$\widehat{\sigma_e} = \mathrm{\bf E}[\theta \, | \, {\bf y}]$', 'FontSize', 28, 'Interpreter', 'latex');
    %end
    xlim([-0.5 8.5]);
    ylim([-0.5 8.5]);
    axis square;
    saveas(f, [dir 'sigmae'], 'epsc');
elseif mode == 6
    % Error vs sigmae with curve for each numiter (ML only)
    if compute
        sigvec = [0.001 0.5:0.5:1 1.75:0.75:10];
        itervec = [1 2 3 50];
        numtrials = 100;
        e_ml_vec = zeros(length(itervec), length(sigvec));
        e_ml_ci = zeros(2*length(itervec), length(sigvec));
        for i = 1 : length(sigvec)
            for k = 1 : length(itervec)
                e_ml = zeros(1, numtrials);
                for j = 1 : numtrials
                    [error_gibbs error_ml] = gibbs_simulation(T, K, N, sigmah, sigvec(i), itervec(k), -1, window, 0);
                    e_ml(j) = error_ml;
                end
                e_ml_vec(k, i) = mean(e_ml);
                e_ml_ci([2*k-1 2*k],i) = bootci(1000, @mean, e_ml);
            end
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    colors = {lighten(lighten('iterml')), lighten('iterml'), 'iterml', [0 0 0]};
    p = [];
    for k = 1 : length(itervec)
        h = plot_errorbar(sigvec, e_ml_vec(k,:), e_ml_ci(2*k-1,:), e_ml_ci(2*k,:), colors{k});
        p = [p h];
    end
    leg = legend(p, '1 iter', '2 iter', '3 iter', '50 iter');  % Make sure to update this
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Noise Level ($\sigma_e)$', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Noise Level', 'FontSize', 14);
    av_1 = mean(e_ml_vec(1,:))
    av_2 = mean(e_ml_vec(2,:))
    av_3 = mean(e_ml_vec(3,:))
    av_50 = mean(e_ml_vec(4,:))
    xlim([-0.5, 11]);
    saveas(f, [dir 'iter_noise'], 'epsc');
elseif mode == 7
    % Error vs N with curve for each numiter (ML only)
    if compute
        Nvec = [10:5:20 30:10:100];
        itervec = [1 2 3 50];
        numtrials = 100;
        e_ml_vec = zeros(length(itervec), length(Nvec));
        e_ml_ci = zeros(2*length(itervec), length(Nvec));
        for i = 1 : length(Nvec)
            for k = 1 : length(itervec)
                e_ml = zeros(1, numtrials);
                for j = 1 : numtrials
                    [error_gibbs error_ml] = gibbs_simulation(T, K, Nvec(i), sigmah, sigmae, itervec(k), -1, window, 0);
                    e_ml(j) = error_ml;
                end
                e_ml_vec(k, i) = mean(e_ml);
                e_ml_ci([2*k-1 2*k],i) = bootci(1000, @mean, e_ml);
            end
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    colors = {lighten(lighten('iterml')), lighten('iterml'), 'iterml', [0 0 0]};
    p = [];
    for k = 1 : length(itervec)
        h = plot_errorbar(Nvec, e_ml_vec(k,:), e_ml_ci(2*k-1,:), e_ml_ci(2*k,:), colors{k});
        p = [p h];
    end
    leg = legend(p, '1 iter', '2 iter', '3 iter', '50 iter');  % Make sure to update this
    set(leg, 'FontSize', 24);
    xlabel('No. Samples ($N$)', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Number of Samples', 'FontSize', 14);
    av_1 = mean(e_ml_vec(1,:))
    av_2 = mean(e_ml_vec(2,:))
    av_3 = mean(e_ml_vec(3,:))
    av_50 = mean(e_ml_vec(4,:))
    xlim([-5 105]);
    saveas(f, [dir 'iter_samples'], 'epsc');
elseif mode == 8
    % f (likelihood) vs numiter
    if compute
        itervec = [2 6 10:10:100];
        numtrials = 1000;
        e_r_vec = zeros(size(itervec));
        e_g_vec = zeros(size(itervec));
        e_ga_vec = zeros(size(itervec));
        e_ml_vec = zeros(size(itervec));
        e_r_ci = zeros(2, length(e_r_vec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(itervec)
            e_r = zeros(1, numtrials);
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                mean_iter = ceil(itervec(i) / 2);  % CHANGE MEAN_ITER HERE
                [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml] = gibbs_simulation(T, K, N, sigmah, sigmae, itervec(i), mean_iter, window, 0);
                e_r(j) = -like_real;
                e_g(j) = -like_g;
                e_ga(j) = -like_ga;
                e_ml(j) = -like_ml;
            end
            e_r_vec(i) = mean(e_r);
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_r_ci(:,i) = bootci(1000, @mean, e_r);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p0 = plot_errorbar(itervec, e_r_vec, e_r_ci(1,:), e_r_ci(2,:), 'real');
    p1 = plot_errorbar(itervec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(itervec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(itervec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p0 p1 p2 p3], 'True Pulses', 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'FontSize', 24, 'Location', 'SouthEast');
    xlabel('No. Iterations ($I$)', 'Interpreter', 'latex', 'FontSize', 16); %24);
    ylabel('$f$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Number of Iterations', 'FontSize', 14);
    xlim([0 105]);
    saveas(f, [dir 'like_iter'], 'epsc');
elseif mode == 9
    % f (likelihood) vs sigmae
    if compute
        sigvec = [0.001 0.5:0.5:6];
        numtrials = 1000;
        e_r_vec = zeros(size(sigvec));
        e_g_vec = zeros(size(sigvec));
        e_ga_vec = zeros(size(sigvec));
        e_ml_vec = zeros(size(sigvec));
        e_r_ci = zeros(2, length(e_r_vec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(sigvec)
            e_r = zeros(1, numtrials);
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2] = gibbs_simulation(T, K, N, sigmah, sigvec(i), numiter, mean_iter, window, 0);
                e_r(j) = -like_real;
                e_g(j) = -like_g;
                e_ga(j) = -like_ga;
                e_ml(j) = -like_ml;
            end
            e_r_vec(i) = mean(e_r);
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_r_ci(:,i) = bootci(1000, @mean, e_r);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p0 = plot_errorbar(sigvec, e_r_vec, e_r_ci(1,:), e_r_ci(2,:), 'real');
    p1 = plot_errorbar(sigvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(sigvec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(sigvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p0 p1 p2 p3], 'True Pulses', 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 16); %24);
    ylabel('$f$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Noise Level', 'FontSize', 14);
    xlim([-0.5, 11]);
    saveas(f, [dir 'like_noise'], 'epsc');
elseif mode == 10
    % Error vs numiter (zoomed in)
    if compute
        itervec = [1 : 15];
        numtrials = 1000;
        %e_g_vec = zeros(size(itervec));
        %e_ga_vec = zeros(size(itervec));
        e_ml_vec = zeros(size(itervec));
        %e_g_ci = zeros(2, length(e_g_vec));
        %e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(itervec)
            disp([int2str(i) ' of ' int2str(length(itervec))]);
            %e_g = zeros(1, numtrials);
            %e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                mean_iter = -1;
                [error_gibbs error_ga error_ml] = gibbs_simulation(T, K, N, sigmah, sigmae, itervec(i), mean_iter, window, 0);
                %e_g(j) = error_gibbs;
                %e_ga(j) = error_ga;
                e_ml(j) = error_ml;
            end
            %e_g_vec(i) = mean(e_g);
            %e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            %e_g_ci(:,i) = bootci(1000, @mean, e_g);
            %e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    %p1 = plot_errorbar(itervec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    %p2 = plot_errorbar(itervec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(itervec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p3], 'IterML');
    set(leg, 'FontSize', 24);
    xlabel('No. Iterations ($I$)', 'Interpreter', 'latex', 'FontSize', 16); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Number of Iterations', 'FontSize', 14);
    xlim([-1 16]);
    saveas(f, [dir 'iter_zoomed'], 'epsc');
elseif mode == 11
    % Error_t vs numiter
    if compute
        itervec = [2 6 10:10:100];
        numtrials = 100;
        e_g_vec = zeros(size(itervec));
        e_ga_vec = zeros(size(itervec));
        e_ml_vec = zeros(size(itervec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(itervec)
            disp([int2str(i) ' of ' int2str(length(itervec))]);
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                mean_iter = ceil(itervec(i) / 2);  % CHANGE MEAN_ITER HERE
                [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2 errort_g errort_ga errort_ml] = gibbs_simulation(T, K, N, sigmah, sigmae, itervec(i), mean_iter, window, 0);
                e_g(j) = errort_g;
                e_ga(j) = errort_ga;
                e_ml(j) = errort_ml;
            end
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);%[prctile(e_g, 5) ; prctile(e_g, 95)];
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(itervec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(itervec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(itervec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'FontSize', 24);
    xlabel('No. Iterations ($I$)', 'Interpreter', 'latex', 'FontSize', 16); %24);
    ylabel('$\mathcal{E}_t$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Number of Iterations', 'FontSize', 14);
    xlim([-3 105]);
    ylim([0 1]);
    saveas(f, [dir 'iter_t'], 'epsc');
elseif mode == 12
    % Error_t vs sigmae
    if compute
        sigvec = [0.001 0.5:0.5:6];
        numtrials = 100;
        e_g_vec = zeros(size(sigvec));
        e_ga_vec = zeros(size(sigvec));
        e_ml_vec = zeros(size(sigvec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        gibbs_all = zeros(0, numtrials);
        for i = 1 : length(sigvec)
            i
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2 errort_g errort_ga errort_ml] = gibbs_simulation(T, K, N, sigmah, sigvec(i), numiter, mean_iter, window, 0);
                e_g(j) = errort_g;
                e_ga(j) = errort_ga;
                e_ml(j) = errort_ml;
            end
            gibbs_all = [gibbs_all ; e_g];
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(sigvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(sigvec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(sigvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 16); %24);
    ylabel('$\mathcal{E}_t$', 'Interpreter', 'latex', 'FontSize', 20);
    %title('Error vs Noise Level', 'FontSize', 14);
    xlim([-0.5, 6.5]);
    ylim([-0.1, 1]);
    saveas(f, [dir 'noise_t'], 'epsc');
elseif mode == 13
    % Error_t vs N (# samples)
    if compute
        Nvec = [10:5:20 30:10:100];
        numtrials = 100;
        e_g_vec = zeros(size(Nvec));
        e_ga_vec = zeros(size(Nvec));
        e_ml_vec = zeros(size(Nvec));
        e_g_ci = zeros(2, length(e_g_vec));
        e_ga_ci = zeros(2, length(e_g_vec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        for i = 1 : length(Nvec)
            e_g = zeros(1, numtrials);
            e_ga = zeros(1, numtrials);
            e_ml = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2 errort_g errort_ga errort_ml] = gibbs_simulation(T, K, Nvec(i), sigmah, sigmae, numiter, mean_iter, window, 0);
                e_g(j) = errort_g;
                e_ga(j) = errort_ga;
                e_ml(j) = errort_ml;
            end
            e_g_vec(i) = mean(e_g);
            e_ga_vec(i) = mean(e_ga);
            e_ml_vec(i) = mean(e_ml);
            e_g_ci(:,i) = bootci(1000, @mean, e_g);
            e_ga_ci(:,i) = bootci(1000, @mean, e_ga);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
        end
        save([dir 'mat/save' int2str(mode)]);
    else
        load([dir 'mat/save' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(Nvec, e_g_vec, e_g_ci(1,:), e_g_ci(2,:), 'gibbs');
    p2 = plot_errorbar(Nvec, e_ga_vec, e_ga_ci(1,:), e_ga_ci(2,:), 'ga');
    p3 = plot_errorbar(Nvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml');
    leg = legend([p1 p2 p3], 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'FontSize', 24);
    xlabel('No. Samples ($N$)', 'Interpreter', 'latex', 'FontSize', 16);% 24);
    ylabel('$\mathcal{E}_t$', 'Interpreter', 'latex', 'FontSize', 20);
    ylim([0 1]);
    xlim([0 105]);
    %title('Error vs Number of Samples', 'FontSize', 14);
    saveas(f, [dir 'samples_t'], 'epsc');
elseif mode == 14
    % Compute runtimes
    e_g_av = 0;
    e_ga_av1 = 0;
    e_ga_av2 = 0;
    e_ml_av = 0;
    e_ann_av = 0;
    numtrials = 100;
    for i = 1 : numtrials
        [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2 errort_g errort_ga errort_ml gibbs_iterTime gaPhase1_iterTime gaPhase2_iterTime ml4_iterTime] = gibbs_simulation(T, K, N, sigmah, sigmae, numiter, mean_iter, window, 0);
        e_g_av = e_g_av + gibbs_iterTime;
        e_ga_av1 = e_ga_av1 + gaPhase1_iterTime;
        e_ga_av2 = e_ga_av2 + gaPhase2_iterTime;
        e_ml_av = e_ml_av + ml4_iterTime;
        ann_time = FRI();  % annihilating filter
        e_ann_av = e_ann_av + ann_time;
    end
    time_g_av = e_g_av / numtrials
    time_ga_av1 = e_ga_av1 / numtrials
    time_ga_av2 = e_ga_av2 / numtrials
    time_ml_av = e_ml_av / numtrials
    time_ann_av = e_ann_av / numtrials
    % Plot bar graphs
    close all;
    f = figure;%('Position', [100 100 1200 400]);
    bar(1:5, [time_ann_av time_g_av time_ga_av1 time_ga_av2 time_ml_av], 0.5, 'FaceColor', 'g', 'LineWidth', 5);
    set(gca, 'XTickLabel', {'Annihilating Filter', 'Gibbs', 'Gibbs-GA', 'Gibbs-GA', 'IterML'});
    ylabel('Runtime for 1 Iteration (sec)', 'FontSize', 16);
    set(gca, 'FontSize', 14);
    saveas(f, [dir 'runtime'], 'epsc');
end

end