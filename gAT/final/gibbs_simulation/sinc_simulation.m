% Run Cadzow vs IterML using periodic sinc kernel

function [error_cad error_ml errort_cad errort_ml t_cad c_cad t_ml c_ml] = sinc_simulation(tau, K, N, sigmae, debug)

% Parameters
%tau = 20;  % [0, tau] is time window containing spikes (20)
%K = 5;  % number of spikes (5)
%N = 30;  % number of samples (30)
%sigmae = 0;  % (2) -- calculated from SNR
%SNR = 5; % (5)
numiter = 3;
window = 10;  % for error calc

if nargin < 5
    debug = 0;
end

T_s = tau/N;
Btau = N;  % odd integer (best when N or N-1)
if mod(N,2) == 0
    Btau = N-1;
end
B = Btau/tau;
M = floor(B*tau/2);

% Random signal
%ck = normrnd(1, 0.4, 1, K);  % old -- not sure why it was like this
ck = normrnd(10, 4, 1, K);
%ck = unifrnd(0, 20, 1, K);
tk = sort(unifrnd(0, tau, 1, K));
delta = 2;
while min(diff(tk)) < delta
    tk = sort(unifrnd(0, tau, 1, K));
end

% Discrepancy here
% Cadzow wants linspace(T_s, tau, N)
% IterML wants linspace(0, tau, N)
% IterML now changed to conform to Cadzow -- see setup_GML_modified
n = linspace(T_s, tau, N);  % points at which samples are taken

% Take samples (sinc kernel)
z = zeros(length(n),1);
for k = 1:length(tk)
    z = z + ck(k)*sinc_kernel(n-tk(k), B, tau)';
end

% Calc sigmae to achieve desired SNR (obsolete -- for testing only)
% TESTING!!!
%start = 1
%tk
%ck
%N
%z
%SNR = 5;
%sigmae = 1;%sqrt(sum(z.^2)/(N*10^(SNR/10))) % now sigmae just taken as input
%SNR = 10*log(sum(z.^2)/(N*sigmae^2))/log(10)

% Add noise
e = normrnd(0, sigmae, length(z), 1);
y = z + e;

%y

% Cadzow reconstruct
[t_cad c_cad time_cad] = run_cadzow_sinc_2(y, tau, K, N, sigmae, tk, ck);



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
[t_ml, c_ml, sigma_ml, elapsed_time, t_inter_ml, c_inter_ml, sig_inter_ml ml4_iterTime] = reconstruct_GML_4(y, tau, args, 0, c_mode, sig_mode);
t_ml = t_ml';
c_ml = c_ml';

% Calculate error
time = linspace(-window, tau+window, 100*N);
zhat_cad = zeros(length(time), 1);
zhat_ml = zeros(length(time), 1);
ztrue = zeros(length(time), 1);
for k = 1 : K
    ztrue = ztrue + ck(k)*sinc_kernel(time-tk(k), B, tau)';
    zhat_cad = zhat_cad + c_cad(k)*sinc_kernel(time-t_cad(k), B, tau)';
    zhat_ml = zhat_ml + c_ml(k)*sinc_kernel(time-t_ml(k), B, tau)';
end
error_cad = (norm(zhat_cad-ztrue)/norm(ztrue))^2;
error_ml = (norm(zhat_ml-ztrue)/norm(ztrue))^2;

errort_cad = t_err(t_cad, tk);
errort_ml = t_err(t_ml, tk);

% Plot results
if debug == 1
    figure;
    title('Noiseless samples (z)');
    plot(z);

    figure;
    title(['Noisy samples (y)']);
    plot(y);

    figure;
    hold on;
    stem(tk,ck,'r');
    stem(t_cad,c_cad,'b');
    stem(t_ml,c_ml,'g');
    legend('True', 'Cadzow', 'IterML');
    title({['error cad = ' num2str(error_cad) ', error ml = ' num2str(error_ml)],['errort cad = ' num2str(errort_cad) ', errort ml = ' num2str(errort_ml)]});
end

% cap error at 1
%error_ml = min(error_ml, 1);
%error_cad = min(error_cad, 1);

end



function e = t_err(t, t_hat)
% Compare spike times

t = sort(t);
t_hat = sort(t_hat);
e = sqrt(sum((t-t_hat).^2)/length(t));

end