% Plot Cramer Rao Bound
% Generate CRB figure for paper 1

function [] = plot_cr()  %type param removed

close all;

compute = 1;  %% compute or load?

% type param (obsolete)
% 0 - gaussian
% 1 - sinc
type = -1;

% mode - string for saving different versions
% '' (empty string) -- original test with sigma 0 to 8
% 'lownoise' -- divided noise by 10
% 'mednoise' -- divided noise by 4
mode = 'mednoise';

T = 20;  % [0, T] is time window containing spikes (20)
K = 5;  % number of spikes (5)
N = 30;  % number of samples (30)
sigmah = 1;
%sigmae = 2;  %2;
numiter = 50;% 50 or 3
mean_iter = 25; % 25 or 2

dir = 'C:\Users\alex\Desktop\Sub-Nyquist Signal Processing\figures\cramer_rao\';

numsigs = 10;
numtrials = 500;

% Run random_signals to generate and save random signals
random_signals(numsigs, K);

if compute == 1

    % choose one
    %if type == 0
        %sigvec = [0.001 0.5:0.5:8];  % for gaussian
    %else
    %    sigvec = [0.01 0.05 0.1 0.2 0.3 0.4 0.6 0.8 1];  % for sinc
    %end
    sigvec = [0.001 1:1:8]/4;  % gaussian sigvec but thinned out
    % TESTING: divided by 10 for smaller noise
    
    % TESTING
    %sigvec = [1];

    e_cr_vec = zeros(1, length(sigvec));
    e_g_vec = zeros(1, length(sigvec));
    e_ga_vec = zeros(1, length(sigvec));
    e_ml_vec = zeros(1, length(sigvec));
    e_cr_psinc_vec = zeros(1, length(sigvec));
    e_cad_vec = zeros(1, length(sigvec));
    e_ml_psinc_vec = zeros(1, length(sigvec));

    t_cr_vec = zeros(1, length(sigvec));
    t_g_vec = zeros(1, length(sigvec));
    t_ga_vec = zeros(1, length(sigvec));
    t_ml_vec = zeros(1, length(sigvec));
    t_cr_psinc_vec = zeros(1, length(sigvec));
    t_cad_vec = zeros(1, length(sigvec));
    t_ml_psinc_vec = zeros(1, length(sigvec));

    for i = 1 : length(sigvec)
        disp([int2str(i) ' of ' int2str(length(sigvec))]);
        sigmae = sigvec(i);
        [devc_cr devc_gibbs devc_GA devc_ml devt_cr devt_gibbs devt_GA devt_ml devc_cr_psinc devt_cr_psinc devc_cad devt_cad devc_ml_psinc devt_ml_psinc] = av_cr(numsigs, numtrials, T, K, N, sigmah, sigmae, numiter, mean_iter, type);
        e_cr_vec(i) = devc_cr;
        e_g_vec(i) = devc_gibbs;
        e_ga_vec(i) = devc_GA;
        e_ml_vec(i) = devc_ml;
        e_cr_psinc_vec(i) = devc_cr_psinc;
        e_cad_vec(i) = devc_cad;
        e_ml_psinc_vec(i) = devc_ml_psinc;
        t_cr_vec(i) = devt_cr;
        t_g_vec(i) = devt_gibbs;
        t_ga_vec(i) = devt_GA;
        t_ml_vec(i) = devt_ml;
        t_cr_psinc_vec(i) = devt_cr_psinc;
        t_cad_vec(i) = devt_cad;
        t_ml_psinc_vec(i) = devt_ml_psinc;
    end
    
    save([dir 'mat/cr_save_' int2str(type) mode]);
    
else
	load([dir 'mat/cr_save_' int2str(type) mode]);
end

%if type == 0

    % Gaussian c
    %{
    f = figure;
    box on;
    hold on;
    p0 = plot_errorbar(sigvec, e_cr_vec, 0, 0, 'real',1);
    p1 = plot_errorbar(sigvec, e_g_vec, 0, 0, 'gibbs',1);
    p2 = plot_errorbar(sigvec, e_ga_vec, 0, 0, 'ga',1);
    p3 = plot_errorbar(sigvec, e_ml_vec, 0, 0, 'iterml',1);
    leg = legend([p0 p1 p2 p3], 'CR', 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('Av Std Dev c', 'FontSize', 24);
    xlim([-0.5, 8.5]);
    ylim([-0.05, 50]);
    saveas(f, [dir 'cr_gauss_c'], 'epsc');
    %}

    % Gaussian t
    f = figure;
    box on;
    hold on;
    p0 = plot_errorbar(sigvec, t_cr_vec, 0, 0, 'real',1);
    p1 = plot_errorbar(sigvec, t_g_vec, 0, 0, 'gibbs',1);
    p2 = plot_errorbar(sigvec, t_ga_vec, 0, 0, 'ga',1);
    p3 = plot_errorbar(sigvec, t_ml_vec, 0, 0, 'iterml',1);
    leg = legend([p0 p1 p2 p3], 'CR', 'Gibbs', 'Gibbs-GA', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('Av Std Dev t', 'FontSize', 24);
    xlim([-0.2, 2.2]);
    ylim([-0.2, 2]);
    saveas(f, [dir 'cr_gauss_t'], 'epsc');

%else

    % Sinc c
    %{
    f = figure;
    box on;
    hold on;
    p4 = plot_errorbar(sigvec, e_cr_psinc_vec, 0, 0, 'real',1);
    p5 = plot_errorbar(sigvec, e_cad_vec, 0, 0, 'cad',1);
    p3 = plot_errorbar(sigvec, e_ml_psinc_vec, 0, 0, 'iterml',1);
    leg = legend([p4 p5 p3], 'CR', 'Cadzow', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('Av Std Dev c', 'FontSize', 24);
    xlim([-0.5, 8.5]);
    ylim([-0.5, 2.5]);
    saveas(f, [dir 'cr_sinc_c'], 'epsc');
    %}

    % Sinc t
    f = figure;
    box on;
    hold on;
    p4 = plot_errorbar(sigvec, t_cr_psinc_vec, 0, 0, 'real',1);
    p5 = plot_errorbar(sigvec, t_cad_vec, 0, 0, 'cad',1);
    p3 = plot_errorbar(sigvec, t_ml_psinc_vec, 0, 0, 'iterml',1);
    leg = legend([p4 p5 p3], 'CR', 'Cadzow', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('Av Std Dev t', 'FontSize', 24);
    xlim([-0.2, 2.2]);
    ylim([-0.2, 2]);
    saveas(f, [dir 'cr_sinc_t'], 'epsc');

%end


end