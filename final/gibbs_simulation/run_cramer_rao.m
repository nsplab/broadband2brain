% Runs Cramer Rao Bound Analysis

function [devc_cr devc_gibbs devc_GA devc_ml devt_cr devt_gibbs devt_GA devt_ml devc_cr_psinc devt_cr_psinc devc_cad devt_cad devc_ml_psinc devt_ml_psinc] = run_cramer_rao(signum, numsigs, numtrials, T, K, N, sigmah, sigmae, numiter, mean_iter, type)

% params (now input)
% type param: see plot_cr
%{
T = 20;  % [0, T] is time window containing spikes
K = 5;  % number of spikes
N = 30;  % number of samples
sigmah = 1;
sigmae = 2;  %2;
numiter = 50;% 50 or 3
mean_iter = 25; % 25 or 2

% trials
numtrials = 10;  % num random trials per true signal (to estimate covariance)
%}

dir = 'C:\Users\alex\Desktop\Sub-Nyquist Signal Processing\figures\cramer_rao\';

% settings
c_mode = 0;  % mode for c update: 0: c_k, 1: c (probably doesn't work anymore)
sig_mode = 0;  % 0: std, 1: E[theta | y]
%delta = 2;  % spike separation

% load random signal
load([dir 'mat/cr_randsigs_' int2str(numsigs)]);
ck = c_mat(signum, :);
tk = t_mat(signum, :);

% results
data_gibbs = zeros(numtrials, 2*K);  % [c | t] vectors
data_GA = zeros(numtrials, 2*K);
data_ml = zeros(numtrials, 2*K);
data_ml_psinc = zeros(numtrials, 2*K);
data_cad = zeros(numtrials, 2*K);

n = linspace(0, T, N);  % points at which samples are taken (gaussian)
n_psinc = linspace(T/N, T, N);  % points at which samples are taken (sinc)

% Take samples (Gaussian kernel)
z = zeros(length(n),1);
for k = 1:length(tk)
    z = z + ck(k)*gausskernel(n-tk(k),sigmah)';
end

% Take samples (periodic sinc kernel)
Btau = N;  % odd integer (best when N or N-1)
if mod(N,2) == 0
    Btau = N-1;
end
tau = T;
B = Btau/tau;
z_cad = zeros(length(n_psinc),1);
for k = 1:length(tk)
    z_cad = z_cad + ck(k)*sinc_kernel(n_psinc-tk(k), B, tau)';
end

% TESTING: print SNR
%start = 1
%tk
%ck
%z_cad
%N
%sigmae
%SNR = 10*log(sum(z_cad.^2)/(N*sigmae^2))/log(10)

M_cr = cramer_rao(n, tk, ck, sigmae, sigmah);
M_cr_psinc = cramer_rao_psinc(n_psinc, tk, ck, sigmae, B, tau);

for i = 1 : numtrials

    % Add noise to samples
    e = normrnd(0, sigmae, length(z), 1);
    y = z + e;
    y_cad = z_cad + e;

    %if type == 0
    
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

        data_gibbs(i, :) = [c_g t_g'];

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

        data_GA(i, :) = [c_ga t_ga'];

        % IterML reconstruct (Gaussian)
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
        [t_ml, c_ml, sigma_ml, time_ml, t_inter_ml, c_inter_ml, sig_inter_ml ml4_iterTime] = reconstruct_GML_4(y, T, args, 0, c_mode, sig_mode);

        data_ml(i, :) = [c_ml' t_ml'];

    %end
    
    %if type == 1
    
        % Cadzow (with sinc kernel)
        [t_cad c_cad time_cad] = run_cadzow_sinc_2(y_cad, tau, K, N, sigmae, tk, ck);

        data_cad(i, :) = [c_cad t_cad];

        % ITERML reconstruct with sinc
        args = struct();
        args.L = N;
        args.K = K;
        args.numiter = numiter;
        args.h = @periodic_sinc;
        args.B = B;
        args.tau = tau;
        args.J = 100;
        c_mode = 0;
        sig_mode = -1;
        %args.delta = 0.2;  % TODO
        method_name = 'GML';
        setup_func = @setup_GML_modified;
        setup_obj = setup_func(tau, 0, args);
        args.setup_obj = setup_obj;
        args.method_name = method_name;
        [t_ml_psinc, c_ml_psinc, sigma_ml, elapsed_time, t_inter_ml, c_inter_ml, sig_inter_ml ml4_iterTime] = reconstruct_GML_4(y_cad, tau, args, 0, c_mode, sig_mode);
        t_ml_psinc = t_ml_psinc';
        c_ml_psinc = c_ml_psinc';

        data_ml_psinc(i, :) = [c_ml_psinc t_ml_psinc];
    
    %end
    
end

cov_gibbs = cov(data_gibbs);
cov_GA = cov(data_GA);
cov_ml = cov(data_ml);
cov_cad = cov(data_cad);
cov_ml_psinc = cov(data_ml_psinc);

% av std dev of c
devc_cr = trace(sqrt(M_cr(1:K,1:K)))/K;
devc_cr_psinc = trace(sqrt(M_cr_psinc(1:K,1:K)))/K;
devc_gibbs = trace(sqrt(cov_gibbs(1:K,1:K)))/K;
devc_GA = trace(sqrt(cov_GA(1:K,1:K)))/K;
devc_ml = trace(sqrt(cov_ml(1:K,1:K)))/K;
devc_cad = trace(sqrt(cov_cad(1:K,1:K)))/K;
devc_ml_psinc = trace(sqrt(cov_ml_psinc(1:K,1:K)))/K;

% av std dev of t
devt_cr = trace(sqrt(M_cr(K+1:2*K,K+1:2*K)))/K;
devt_cr_psinc = trace(sqrt(M_cr_psinc(K+1:2*K,K+1:2*K)))/K;
devt_gibbs = trace(sqrt(cov_gibbs(K+1:2*K,K+1:2*K)))/K;
devt_GA = trace(sqrt(cov_GA(K+1:2*K,K+1:2*K)))/K;
devt_ml = trace(sqrt(cov_ml(K+1:2*K,K+1:2*K)))/K;
devt_cad = trace(sqrt(cov_cad(K+1:2*K,K+1:2*K)))/K;
devt_ml_psinc = trace(sqrt(cov_ml_psinc(K+1:2*K,K+1:2*K)))/K;

% better way?: take std with true value as "mean"

%devt_cr_psinc
%devt_cad

% TESTING
%{
figure;
hold on;
stem(tk, ck, 'c');
stem(t_cad, c_cad, 'k');
stem(t_ml, c_ml, 'g');
%}
