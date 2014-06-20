% Runs Gibbs sampling and iterative ML
% Single run with parameters taken as inputs
% Generates plots if debug set to 1
% Set mean_iter = -1 to not run gibbs

function [error_g error_ga error_ml sigma_g sigma_ml like_real like_g like_ga like_ml flag1 flag2 errort_g errort_ga errort_ml gibbs_iterTime gaPhase1_iterTime gaPhase2_iterTime ml4_iterTime] = gibbs_simulation(T, K, N, sigmah, sigmae, numiter, mean_iter, window, debug, sig_mode)

% settings
c_mode = 0;  % 0: c_k, 1: c (probably doesn't work anymore)
if nargin < 10
    sig_mode = 0;  % 0: std, 1: E[theta | y]
end
spike_mode = 1;  % 0: single example, 1: random, 2: fig1, 3: fig9
delta = 2;  % spike separation

% Generate spikes
if spike_mode == 0
    ck = [6 12 10 15 9];%   13 16 9 10 11];  % hard-coded spikes used by Tan & Goyal
    tk = [9 13 16 3 6];%   14 18 7 10 5];
    K = length(ck);
elseif spike_mode == 1
    ck = normrnd(10, 4, 1, K);
    %ck = unifrnd(0, 20, 1, K);
    tk = sort(unifrnd(0, 20, 1, K));
    while min(diff(tk)) < delta
        tk = sort(unifrnd(0, 20, 1, K));
    end
elseif spike_mode == 2  % figure 1: flowchart
    K = 2;
    tk = [2 5];
    ck = [7 6];
    T = 8;
    N = 6;
    sigmae = 2;
elseif spike_mode == 3  % figure 9: Gaussian kernel fail
    tk = [5.6892 12.8557 14.0082 18.0777 19.8478]-3;
    ck = [12.2782 11.0373 9.3133 6.1440 4.6548];
elseif spike_mode == 4  % negative spikes
    ck = normrnd(0, 10, 1, K);
    tk = sort(unifrnd(0, 20, 1, K));
    while min(diff(tk)) < delta
        tk = sort(unifrnd(0, 20, 1, K));
    end
elseif spike_mode == 5  % copy genetic algorithm paper
    T = 14;
    t_prev = 0;
    tk= [];
    while t_prev <= 14
        t_prev = t_prev + 0.5 + exprnd(1);
        if t_prev <= 14
            tk = [tk t_prev];
        end
    end
    K = length(tk);
    ck = unifrnd(3, 10, 1, K);
    %N = 100;
    %sigmae = sqrt(0.25);
    sigmah = T / (N-1);
    % 1000 iterations
end

n = linspace(0, T, N);  % points at which samples are taken

% Take samples (Gaussian kernel)
z = zeros(length(n),1);
for k = 1:length(tk)
    z = z + ck(k)*gausskernel(n-tk(k),sigmah)';
end

inc = -1;

% loop until satisfactory E
%inc = 0.025;
%target_E = 0.04;
%tol = 0.001;
%error_ml = -1;
%error_g = -1;
%error_ga = -1;
%while abs(error_ga - target_E) > tol
    
    %error_ml

    % Add noise to samples
    e = normrnd(0, sigmae, length(z), 1);
    y = z + e;

    % Gibbs reconstruct
    args = struct();
    args.sigmah = sigmah;
    args.L = N;
    args.K = K;
    args.numiter = numiter;
    args.mean_iter = mean_iter;
    if mean_iter == -1
        % hack to fix bug
        % don't need results, just need this to not throw an error
        % for running ML with numiter = 1 (actually running with mean_iter
        % = 1 works now)
        args.numiter = 2;
        args.mean_iter = 1;
    end
    [t_g, c_g, sigma_g, elapsed_time, t_inter_g, c_inter_g, sig_inter_g gibbs_iterTime] = reconstruct_gibbs(y, T, args);
    %gibbs_time = elapsed_time
    args_g = args;

    % GA reconstruct
    args = struct();
    args.sigmah = sigmah;
    args.L = N;
    args.K = K;
    args.I_e = floor(numiter/2);  % TESTING: switching these
    args.I_m = ceil(numiter/2);
    args.sigmae = sigmae;
    args.t = tk;
    args.c = ck;
    [t_ga, c_ga, sigma_ga, elapsed_time_ga, t_inter_ga, c_inter_ga, sig_inter_ga gaPhase1_iterTime gaPhase2_iterTime] = reconstruct_GA(y, T, args);

    % IterML reconstruct
    args = struct();
    args.L = N;
    args.K = K;
    args.numiter = numiter;
    args.h = @gaussian;
    args.sigmah = sigmah;
    args.J = 100;
    %args.delta = 0.2;  % TODO
    method_name = 'GML';
    setup_func = eval(['@setup_' method_name]);
    setup_obj = setup_func(T, 0, args);
    args.setup_obj = setup_obj;
    args.method_name = method_name;
    [t_ml, c_ml, sigma_ml, elapsed_time, t_inter_ml, c_inter_ml, sig_inter_ml ml4_iterTime] = reconstruct_GML_4(y, T, args, 0, c_mode, sig_mode);  %TODO: change GML variant
    %[t_ml, c_ml, sigma_ml, elapsed_time, t_inter_ml, c_inter_ml, sig_inter_ml] = reconstruct_GML_6(y, T, args, 0, c_mode, sig_mode);  %TODO: change GML variant
    %ml4_iterTime = 0;
    %[t_ml5, c_ml5, sigma_ml5, elapsed_time5, t_inter_ml5, c_inter_ml5, sig_inter_ml5] = reconstruct_GML_5(y, T, args, 0, c_mode, sig_mode);  %TODO: change GML variant
    %[t_ml6, c_ml6, sigma_ml6, elapsed_time6, t_inter_ml6, c_inter_ml6, sig_inter_ml6] = reconstruct_GML_6(y, T, args, 0, c_mode, sig_mode);  %TODO: change GML variant
    ml_time = elapsed_time;

    % Calculate error
    time = linspace(-window, T+window, 100*N);
    zhat_g = zeros(length(time), 1);
    zhat_ga = zeros(length(time), 1);
    zhat_ml = zeros(length(time), 1);
    %zhat_ml5 = zeros(length(time), 1);
    %zhat_ml6 = zeros(length(time), 1);
    ztrue = zeros(length(time), 1);
    for k = 1 : K
        ztrue = ztrue + ck(k)*gausskernel(time-tk(k), sigmah)';
        zhat_g = zhat_g + c_g(k)*gausskernel(time-t_g(k), sigmah)';
        zhat_ga = zhat_ga + c_ga(k)*gausskernel(time-t_ga(k), sigmah)';
        zhat_ml = zhat_ml + c_ml(k)*gausskernel(time-t_ml(k), sigmah)';
    %    zhat_ml5 = zhat_ml5 + c_ml5(k)*gausskernel(time-t_ml5(k), sigmah)';
    %    zhat_ml6 = zhat_ml6 + c_ml6(k)*gausskernel(time-t_ml6(k), sigmah)';
    end
    error_g = (norm(zhat_g-ztrue)/norm(ztrue))^2;
    error_ga = (norm(zhat_ga-ztrue)/norm(ztrue))^2;
    error_ml = (norm(zhat_ml-ztrue)/norm(ztrue))^2;
    %error_ml5 = (norm(zhat_ml5-ztrue)/norm(ztrue))^2;
    %error_ml6 = (norm(zhat_ml6-ztrue)/norm(ztrue))^2;
    
    errort_g = t_err(t_g, tk);
    errort_ga = t_err(t_ga, tk);
    errort_ml = t_err(t_ml, tk);

%end  % for loop until target_E

%error_g
%error_ga
%error_ml
    
% Test whether ML zhat or y approximates z better
%{
zhat_ml_s = zeros(length(n),1);
for k = 1:length(tk)
    zhat_ml_s = zhat_ml_s + c_ml(k)*gausskernel(n-t_ml(k),sigmah)';
end
error_ml_s = (norm(zhat_ml_s-z)/norm(z))^2
error_y_s = (norm(y-z)/norm(z))^2
%}

% Test whether ML achieves better likelihood than real spikes
%if debug
    
    zhat_ml_s = zeros(length(n),1);
    for k = 1:length(tk)
        zhat_ml_s = zhat_ml_s + c_ml(k)*gausskernel(n-t_ml(k),sigmah)';
    end
%{
    zhat_ml_s5 = zeros(length(n),1);
    for k = 1:length(tk)
        zhat_ml_s5 = zhat_ml_s5 + c_ml5(k)*gausskernel(n-t_ml5(k),sigmah)';
    end
    
    zhat_ml_s6 = zeros(length(n),1);
    for k = 1:length(tk)
        zhat_ml_s6 = zhat_ml_s6 + c_ml6(k)*gausskernel(n-t_ml6(k),sigmah)';
    end
%}    
    zhat_g_s = zeros(length(n),1);
    for k = 1:length(tk)
        zhat_g_s = zhat_g_s + c_g(k)*gausskernel(n-t_g(k),sigmah)';
    end
    
    zhat_ga_s = zeros(length(n),1);
    for k = 1:length(tk)
        zhat_ga_s = zhat_ga_s + c_ga(k)*gausskernel(n-t_ga(k),sigmah)';
    end
    
    like_real = sum((z-y).^2);
    like_g = sum((zhat_g_s-y).^2);
    like_ga = sum((zhat_ga_s-y).^2);
    like_ml = sum((zhat_ml_s-y).^2);
    %like_ml5 = sum((zhat_ml_s5-y).^2)
    flag1 = 1;
    flag2 = 1;
%end

if debug == 1
    
    if spike_mode == 2
        
        dir = 'C:\Users\alex\Desktop\figures\fig1\';
        
        % x(t)
        f = fig1settings;
        plot_stem(tk, ck, 'real', 0, 1);
        saveas(f, [dir 'x'], 'epsc');
        
        % z(t)
        f = fig1settings;
        plot_stem(tk, ck, lighten('real',1), 0, 1);
        plot_line(time, ztrue, 'real', 0, 1);
        saveas(f, [dir 'z'], 'epsc');
        
        % z[n]
        f = fig1settings;
        plot_line(time, ztrue, lighten('real',1), 0, 1);
        plot_dots(n, z, 'real', 1);
        saveas(f, [dir 'zn'], 'epsc');
        
        % y[n]
        f = fig1settings;
        plot_line(time, ztrue, lighten('real',1), 0, 1);
        plot_stem(n, y, 'err', 0, 1, z);
        %for i = 1 : length(n)
        %    plot_line([n(i) n(i)], [z(i) y(i)], 'err');
        %end
        plot_dots(n, z, lighten('real',1), 1);
        plot_dots(n, y, 'err', 1);
        saveas(f, [dir 'y'], 'epsc');
        
        % x_hat(t)
        f = fig1settings;
        plot_stem(tk, ck, lighten('real',1), 0, 1);
        plot_stem(t_g, c_g, 'gibbs', 0, 1);
        saveas(f, [dir 'xhat'], 'epsc');
        
        % z_hat(t)
        f = fig1settings;
        plot_stem(t_g, c_g, lighten('gibbs',1), 0, 1);
        plot_line(time, ztrue, lighten('real',1), 0, 1);
        plot_line(time, zhat_g, 'gibbs', 0, 1);
        saveas(f, [dir 'zhat'], 'epsc');
        
        % E
        f = fig1settings;
        jbfill(time, max(ztrue, zhat_g)', min(ztrue, zhat_g)', lighten('err',1));
        hold on;
        plot_line(time, ztrue, 'real', 0, 1);
        plot_line(time, zhat_g, 'gibbs', 0, 1);
        saveas(f, [dir 'E'], 'epsc');
        
        % Text
        f = figure;
        title(['$x(t) z(t) h(t) z[n] y[n] e[n] \widehat{x}(t) \widehat{z}(t) \sigma_e \mathcal{E} = $' num2str(error_g)], 'FontSize', 18, 'Interpreter', 'latex');
        set(gca, 'FontSize', 14);
        xlim([0 5]);
        ylim([0 0.5]);
        saveas(f, [dir 'text'], 'epsc');
        
        f = figure;
        set(gca, 'FontSize', 14);
        xlim([0 0.1]);
        set(gca, 'XTick', 0 : 0.025 : 0.1);
        saveas(f, [dir 'text2'], 'epsc');
        
        f = figure;
        set(gca, 'FontSize', 14);
        xlim([0 0.5]);
        set(gca, 'XTick', 0 : 0.025 : 0.5);
        saveas(f, [dir 'text3'], 'epsc');
        
    else
        
        dir = 'C:\Users\alex\Desktop\figures\';
        
        % x
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 36);
        plot_stem(tk, ck, 'real');
        plot([-5 25], [0 0], 'k');
        axis([-5 25 -5 20]);
        xlabel('Time', 'FontSize', 45);
        ylabel('Amplitude', 'FontSize', 45);
        saveas(f, [dir 'x'], 'epsc');
        
        % Gibbs z
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, 'real');
        plot_line(time, zhat_g, 'gibbs');
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        title(['Gibbs: $\varepsilon = ' num2str(error_g) '$'], 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend('${\bf z(t)}$', '${\bf \widehat{z}(t)}$');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        %saveas(f, [dir 'gibbs_z'], 'epsc');

        % Gibbs x
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        p1 = plot_stem(tk, ck, 'real');
        p2 = plot_stem(t_g, c_g, 'gibbs');
        plot([-5 25], [0 0], 'k');
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        title('Gibbs', 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend([p1 p2], '${\bf x(t)}$', '${\bf \widehat{x}(t)}$');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        %saveas(f, [dir 'gibbs_x'], 'epsc');
    
        % Gibbs-GA z
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, 'real');
        plot_line(time, zhat_ga, 'ga');
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        %title(['Gibbs: $\varepsilon = ' num2str(error_g) '$'], 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend('${\bf z(t)}$', '${\bf \widehat{z}(t)}$');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'gibbs_z'], 'epsc');

        % Gibbs-GA x
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        p1 = plot_stem(tk, ck, 'real');
        p2 = plot_stem(t_ga, c_ga, 'ga');
        plot([-5 25], [0 0], 'k');
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        %title('Gibbs', 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend([p1 p2], '${\bf x(t)}$', '${\bf \widehat{x}(t)}$');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'gibbs_x'], 'epsc');

        % ML z
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, 'real');
        plot_line(time, zhat_ga, 'ga');
        plot_line(time, zhat_ml, 'iterml');
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        %title(['IterML: $\varepsilon = ' num2str(error_g) '$'], 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend('${\bf z(t)}$', '${\bf \widehat{z}(t)}$ (Gibbs-GA)', '${\bf \widehat{z}(t)}$ (IterML)');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'iterml_z'], 'epsc');

        % ML x
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        p1 = plot_stem(tk, ck, 'real');
        p2 = plot_stem(t_ga, c_ga, 'ga');
        p3 = plot_stem(t_ml, c_ml, 'iterml');
        %plot(n, y, 'm.-');
        plot([-5 25], [0 0], 'k');
        axis([-5 25 -5 20]);
        %title('IterML', 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend([p1 p2 p3], '${\bf x(t)}$', '${\bf \widehat{x}(t)}$ (Gibbs-GA)', '${\bf \widehat{x}(t)}$ (IterML)');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'iterml_x'], 'epsc');
        
        % ML (only) z
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, 'real');
        %plot_line(time, zhat_g, lighten('gibbs',1));
        plot_line(time, zhat_ml, 'iterml', 1);
        %plot(n, y, 'm.-');
        axis([-5 25 -5 20]);
        %title(['IterML: $\varepsilon = ' num2str(error_g) '$'], 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend('${\bf z(t)}$', '${\bf \widehat{z}(t)}$ (IterML)');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'mlonly_z'], 'epsc');

        % ML (only) x
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        p1 = plot_stem(tk, ck, 'real');
        %p2 = plot_stem(t_g, c_g, lighten('gibbs',1));
        p3 = plot_stem(t_ml, c_ml, 'iterml', 1);
        %plot(n, y, 'm.-');
        plot([-5 25], [0 0], 'k');
        axis([-5 25 -5 20]);
        %title('IterML', 'FontSize', 14, 'Interpreter', 'latex');
        leg = legend([p1 p3], '${\bf x(t)}$', '${\bf \widehat{x}(t)}$ (IterML)');
        set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 28);
        ylabel('Amplitude', 'FontSize', 28);
        saveas(f, [dir 'mlonly_x'], 'epsc');
        
        % E (error snippet)
        f = figure;
        jbfill(time, max(ztrue, zhat_ml)', min(ztrue, zhat_ml)', lighten('err',1));
        hold on;
        box on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, 'real');
        plot_line(time, zhat_ml, 'iterml');
        axis([-5 25 -5 20]);
        set(gca, 'XTick', [], 'YTick', []);
        if inc ~= -1
            saveas(f, [dir '/fig5/E' num2str(inc*round(error_ml/inc))], 'epsc');
        end
        
        % x snippet
        f = figure;
        box on;
        hold on;
        set(gca, 'FontSize', 18);
        p1 = plot_stem(tk, ck, 'real');
        %p2 = plot_stem(t_g, c_g, lighten('gibbs',1));
        p3 = plot_stem(t_ml, c_ml, 'iterml');
        %plot(n, y, 'm.-');
        plot([-5 25], [0 0], 'k');
        axis([-5 25 -5 20]);
        set(gca, 'XTick', [], 'YTick', []);
        saveas(f, [dir '/fig5/x' num2str(inc*round(error_ml/inc))], 'epsc');
        
        % y[n] (noise snippet)
        f = figure;
        hold on;
        box on;
        set(gca, 'FontSize', 18);
        plot_line(time, ztrue, lighten('real',1));
        plot_stem(n, y, 'err', 0, 0, z);
        plot_dots(n, z, 'real');
        plot_dots(n, y, 'err');
        axis([-5 25 -5 20]);
        set(gca, 'XTick', [], 'YTick', []);
        saveas(f, [dir '/fig5/y-sig' num2str(sigmae)], 'epsc');
        
        %{
        % ML 5 z
        f = figure;
        box on;
        hold on;
        plot_line(time, ztrue, 'real');
        plot_line(time, zhat_ml, 'iterml', 1);
        plot_line(time, zhat_ml5, 'err', 1);
        plot(n, y, 'm.-');
        axis([-5 25 -20 20]);
        %title(['IterML: $\varepsilon = ' num2str(error_g) '$'], 'FontSize', 14, 'Interpreter', 'latex');
        %leg = legend('${\bf z(t)}$', '4', '5');
        %set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 24);
        ylabel('Amplitude', 'FontSize', 24);
        saveas(f, [dir 'mlonly_z'], 'epsc');

        % ML 5 x
        f = figure;
        box on;
        hold on;
        p1 = plot_stem(tk, ck, 'real');
        p2 = plot_stem(t_ml, c_ml, 'iterml', 1);
        p3 = plot_stem(t_ml5, c_ml5, 'err', 1);
        %plot(n, y, 'm.-');
        axis([-5 25 -20 20]);
        %title('IterML', 'FontSize', 14, 'Interpreter', 'latex');
        %leg = legend([p1 p2 p3], '${\bf x(t)}$', '4', '5');
        %set(leg, 'Interpreter', 'latex', 'FontSize', 24);
        xlabel('Time', 'FontSize', 24);
        ylabel('Amplitude', 'FontSize', 24);
        saveas(f, [dir 'mlonly_x'], 'epsc');
        %}

        % Plot convergence of params

        %{
        % Gibbs t
        figure;
        box on;
        hold on;
        for k = 1 : K
            h = plot([1 numiter*1.2], [tk(k) tk(k)], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            set(h, 'LineWidth', 2);
        end
        for k = 1 : K
            h = plot(numiter*1.1, t_g(k), 'dk');
            set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        end
        h = plot((numiter-mean_iter+0.5)*ones(1,2), [-5, 20], 'k--');
        set(h, 'LineWidth', 2);
        plot(t_inter_g', 'LineWidth', 2);
        axis([1 numiter*1.2 -5 20]);
        title('Gibbs ${\bf t}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf t}$', 'FontSize', 12, 'Interpreter', 'latex');

        % Gibbs c
        figure;
        box on;
        hold on;
        for k = 1 : K
            h = plot([1 numiter*1.2], [ck(k) ck(k)], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            set(h, 'LineWidth', 2);
        end
        for k = 1 : K
            h = plot(numiter*1.1, c_g(k), 'dk');
            set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        end
        h = plot((numiter-mean_iter+0.5)*ones(1,2), [-5, 20], 'k--');
        set(h, 'LineWidth', 2);
        plot(c_inter_g', 'LineWidth', 2);
        axis([1 numiter*1.2 -5 20]);
        title('Gibbs ${\bf c}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf c}$', 'FontSize', 12, 'Interpreter', 'latex');

        % Gibbs sigmae
        figure;
        box on;
        hold on;
        h = plot([1 numiter*1.2], [sigmae sigmae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        set(h, 'LineWidth', 2);
        h = plot((numiter-mean_iter+0.5)*ones(1,2), [-5, 20], 'k--');
        set(h, 'LineWidth', 2);
        h = plot(sig_inter_g');
        set(h, 'LineWidth', 2);
        h = plot(numiter*1.1, sigma_g, 'dk');
        set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        axis([1 numiter*1.2 -1 2*sigmae+1]);
        title('Gibbs ${\bf \sigma_e}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf \sigma_e}$', 'FontSize', 12, 'Interpreter', 'latex');

        % IterML t
        figure;
        box on;
        hold on;
        for k = 1 : K
            h = plot([1 numiter*1.2], [tk(k) tk(k)], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            set(h, 'LineWidth', 2);
        end
        for k = 1 : K
            h = plot(numiter*1.1, t_ml(k), 'dk');
            set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        end
        plot(t_inter_ml', 'LineWidth', 2);
        axis([1 numiter*1.2 -5 20]);
        title('IterML ${\bf t}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf t}$', 'FontSize', 12, 'Interpreter', 'latex');

        % IterML c
        figure;
        box on;
        hold on;
        for k = 1 : K
            h = plot([1 numiter*1.2], [ck(k) ck(k)], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            set(h, 'LineWidth', 2);
        end
        for k = 1 : K
            h = plot(numiter*1.1, c_ml(k), 'dk');
            set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        end
        plot(c_inter_ml', 'LineWidth', 2);
        axis([1 numiter*1.2 -5 20]);
        title('IterML ${\bf c}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf c}$', 'FontSize', 12, 'Interpreter', 'latex');

        % IterML sigmae
        figure;
        box on;
        hold on;
        h = plot([1 numiter*1.2], [sigmae sigmae], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        set(h, 'LineWidth', 2);
        %h = plot(sig_inter_ml');
        %set(h, 'LineWidth', 2);
        h = plot(numiter*1.1, sigma_ml, 'dk');
        set(h, 'MarkerSize', 5, 'MarkerFaceColor', 'k');
        axis([1 numiter*1.2 -1 2*sigmae+1]);
        title('IterML ${\bf \sigma_e}$', 'FontSize', 14, 'Interpreter', 'latex');
        xlabel('Iteration', 'FontSize', 12);
        ylabel('${\bf \sigma_e}$', 'FontSize', 12, 'Interpreter', 'latex');
        %}
        
    end
    
end

% cap error at 1
error_g = min(error_g, 1);
error_ga = min(error_ga, 1);
error_ml = min(error_ml, 1);
    
end

function e = t_err(t, t_hat)
% Compare spike times

t = sort(t);
t_hat = sort(t_hat');
e = sqrt(sum((t-t_hat).^2)/length(t));

end