% Tests fast nonnegative deconvolution and iter ML on simulated calcium
% imaging data

function [runtime_dec runtime_ml] = test_cal(T, lambda)

est_params = 0;

% Parameters (for simulated data)
if nargin < 1
    T = 400;
end
delta = 0.0333;
alpha = 1;
beta = 0; %0.5;
sigma = 0.2;
tau = 1;
if nargin < 2
    lambda = 1;
end
gamma = 1 - delta*tau;
time = delta : delta : T*delta;

V.T = T;
V.dt = delta;

% Simulate data
[n C F] = sim_cal_data(T, delta, lambda, alpha, beta, sigma, gamma);

% Run deconvolution
if est_params
    V.est_sig = 1;
    V.est_lam = 1;
    V.est_gam = 1;
    V.est_b = 1;
    V.est_a = 1;
    V.fast_iter_max = 50;  % why does this break if > 1000
    P = [];
else
    P.a = alpha;
    P.b = beta;
    P.sig = sigma;
    P.gam = 1-delta/tau;
    P.lam = lambda*delta;
end
ticID = tic;
[n_dec P_dec V_dec C_dec] = fast_oopsi(F, V, P);
runtime_dec = toc(ticID);

% Run iter ML
args = struct();
args.L = T;
args.K = ceil(1.2*lambda*delta*T);  % TODO
args.K
args.numiter = 7;
args.h = @exponential;
args.tau = tau;
args.a = alpha;
args.t_resolution = delta;  % TODO
args.delta = delta;  % TODO
args.delta_t = get_dt();
method_name = 'RMSE_cal';
args.method_name = method_name;
setup_func = eval(['@setup_' method_name]);
setup_obj = setup_func(T*delta, 0, args);
args.setup_obj = setup_obj;
args.method_name = method_name;
[t_ml, c_ml, sigma_ml, elapsed_time] = reconstruct_RMSE_cal(F', T*delta, args);
%[t_ml, c_ml, sigma_ml, elapsed_time] = reconstruct_GML_4(F', T*delta, args);
C_ml = model_samples(0, time, T*delta, args, [], t_ml, c_ml, sigma_ml);
runtime_ml = elapsed_time;

% Plot results
close all;

% Spikes
figure;
hold on;
stem(time, n, 'c');
stem(time, n_dec*P_dec.a+P_dec.b, 'b--');
stem(t_ml, c_ml, 'k');
title('n');
legend('real', 'deconvolution', 'ML');

% Calcium
figure;
hold on;
plot(time, C, 'c');
plot(time, F, 'r--');
plot(time, C_dec*P_dec.a+P_dec.b, 'b--');
plot(time, C_ml, 'k--');
title('C');
legend('real', 'noisy', 'deconvolution', 'ML');

end