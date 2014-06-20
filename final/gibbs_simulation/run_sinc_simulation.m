% Run Cadzow vs IterML using periodic sinc kernel
% Generate figure for paper 1
% Note: see cadzow_denoise for iteration settings

function [] = run_sinc_simulation(mode, compute)

close all;

% Parameters
tau = 20;  % [0, tau] is time window containing spikes (20)
K = 5;  % number of spikes (5)
N = 30;  % number of samples (30)
sigmae = 0.8;  % (0.4)
numiter = 3;  % (3)
window = 10;  % for error calc

%sigvec = [0.01 0.05 0.1 0.2 0.3 0.4 0.6 0.8 1];  % for mode = 3 or 5
sigvec = [0.001 0.5:0.5:8];  % this is the gaussian sigvec

% Note that which type of error to use (E or t) is hard-coded later

dir = 'C:\Users\alex\Desktop\Sub-Nyquist Signal Processing\figures\';

% CHANGE MODE HERE
if nargin < 1
    mode = 3;
    compute = 1;  % 1 to compute, 0 to load from .mat files
end

if mode == 0
    % Single trial
    [error_cad error_ml errort_cad errort_ml] = sinc_simulation(tau, K, N, sigmae, 1);
elseif mode == 3
    % Error vs sigmae
    if compute
        numtrials = 10;
        e_ml_vec = zeros(size(sigvec));
        e_cad_vec = zeros(size(sigvec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        e_cad_ci = zeros(2, length(e_cad_vec));
        for i = 1 : length(sigvec)
            i
            e_ml = zeros(1, numtrials);
            e_cad = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_cad error_ml errort_cad errort_ml] = sinc_simulation(tau, K, N, sigvec(i));
                % Decide which type of error to use
                e_ml(j) = error_ml;
                e_cad(j) = error_cad;
            end
            e_ml_vec(i) = mean(e_ml);
            e_cad_vec(i) = mean(e_cad);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
            e_cad_ci(:,i) = bootci(1000, @mean, e_cad);
        end
        save([dir 'mat/save_sinc' int2str(mode)]);
    else
        load([dir 'mat/save_sinc' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(sigvec, e_cad_vec, e_cad_ci(1,:), e_cad_ci(2,:), 'cad',1);  % last param to 1 for no errorbars
    p2 = plot_errorbar(sigvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml',1);
    leg = legend([p1 p2], 'Cadzow', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    %title('Error vs Noise Level', 'FontSize', 14);
    xlim([-0.5, 8.5]);
    ylim([-0.3, 2.7]);
    saveas(f, [dir 'sinc_noise'], 'epsc');
    
elseif mode == 4
    % error vs num samples
    if compute
        Nvec = [10:5:20 30:10:100];
        Nvec(1) = 11;  % for Cadzow
        numtrials = 100;
        e_ml_vec = zeros(size(Nvec));
        e_cad_vec = zeros(size(Nvec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        e_cad_ci = zeros(2, length(e_cad_vec));
        for i = 1 : length(Nvec)
            i
            e_ml = zeros(1, numtrials);
            e_cad = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_cad error_ml errort_cad errort_ml] = sinc_simulation(tau, K, Nvec(i), sigmae);
                % Decide which type of error to use (E or t)
                e_ml(j) = errort_ml;
                e_cad(j) = errort_cad;
            end
            e_ml_vec(i) = mean(e_ml);
            e_cad_vec(i) = mean(e_cad);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
            e_cad_ci(:,i) = bootci(1000, @mean, e_cad);
        end
        save([dir 'mat/save_sinc' int2str(mode)]);
    else
        load([dir 'mat/save_sinc' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(Nvec, e_cad_vec, e_cad_ci(1,:), e_cad_ci(2,:), 'cad',1);
    p2 = plot_errorbar(Nvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml',1);
    leg = legend([p1 p2], 'Cadzow', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('No. Samples ($N$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}$', 'Interpreter', 'latex', 'FontSize', 32);
    %title('Error vs Noise Level', 'FontSize', 14);
    %xlim([-0.5, 8.5]);
    %ylim([-0.05, 0.8]);
    saveas(f, [dir 'sinc_noise'], 'epsc');
    
elseif mode == 5
    % ErrorT vs sigmae
    if compute
        numtrials = 1000;
        e_ml_vec = zeros(size(sigvec));
        e_cad_vec = zeros(size(sigvec));
        e_ml_ci = zeros(2, length(e_ml_vec));
        e_cad_ci = zeros(2, length(e_cad_vec));
        for i = 1 : length(sigvec)
            i
            e_ml = zeros(1, numtrials);
            e_cad = zeros(1, numtrials);
            for j = 1 : numtrials
                [error_cad error_ml errort_cad errort_ml] = sinc_simulation(tau, K, N, sigvec(i));
                % Decide which type of error to use
                e_ml(j) = errort_ml;
                e_cad(j) = errort_cad;
            end
            e_ml_vec(i) = mean(e_ml);
            e_cad_vec(i) = mean(e_cad);
            e_ml_ci(:,i) = bootci(1000, @mean, e_ml);
            e_cad_ci(:,i) = bootci(1000, @mean, e_cad);
        end
        save([dir 'mat/save_sinc' int2str(mode)]);
    else
        load([dir 'mat/save_sinc' int2str(mode)]);
    end
    f = figure;
    box on;
    hold on;
    p1 = plot_errorbar(sigvec, e_cad_vec, e_cad_ci(1,:), e_cad_ci(2,:), 'cad',0);  % last param to 1 for no errorbars
    p2 = plot_errorbar(sigvec, e_ml_vec, e_ml_ci(1,:), e_ml_ci(2,:), 'iterml',0);
    leg = legend([p1 p2], 'Cadzow', 'IterML');
    set(leg, 'Location', 'NorthWest', 'FontSize', 24);
    set(gca, 'FontSize', 18);
    xlabel('Noise Level ($\sigma_e$)', 'Interpreter', 'latex', 'FontSize', 24); %24);
    ylabel('$\mathcal{E}_t$', 'Interpreter', 'latex', 'FontSize', 32);
    %title('Error vs Noise Level', 'FontSize', 14);
    xlim([-0.5, 8.5]);
    ylim([-0.3, 3.3]);
    saveas(f, [dir 'sinc_noise_t'], 'epsc');
    
    
end




end