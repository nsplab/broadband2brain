% Implementation of Blu "Sparse Sampling of Signal Innovations" with sinc
% kernel
% Gaussian kernel seeme unstable so I'm trying to reproduce Blu's results
% exactly

close all;
clear all;

% Parameters
tau = 20;  % [0, tau] is time window containing spikes (20)
K = 5;  % number of spikes (5)
N = 30;  % number of samples (30)
%sigmae = 0;  % (2)
SNR = 5; % (5)
numiter = 3;

T_s = tau/N;
Btau = N;  % odd integer (best when N or N-1)
if mod(N,2) == 0
    Btau = N-1;
end
B = Btau/tau;
M = floor(B*tau/2);

B
tau

% Random signal
ck = normrnd(1, 0.4, 1, K);
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

% Calc sigmae to achieve desired SNR
sigmae = sqrt(sum(z.^2)/(N*10^(SNR/10)))

% Add noise
e = normrnd(0, sigmae, length(z), 1);
y = z + e;

%y

%%% Beginning of reconstruction algorithm

ticID = tic;

% Compute DFT y_hat (can apparently be made faster using FFT)
orig_ID = tic;
y_hat = zeros(1, Btau);
clear i;
for m = -M : M
    for nn = 1 : N
        y_hat(m+M+1) = y_hat(m+M+1) + y(nn)*exp(-i*2*pi*nn*m/N);
    end
end
orig_time = toc(orig_ID)

fft_ID = tic;
y_fft = fft(y, Btau)';
fft_time = toc(fft_ID)

%y_hat - y_fft


%y_hat

% Cadzow denoising step
y_hat = cadzow_denoise(y_hat, K);

%y_hat

% Form matrix S
S = zeros(2*M-K, K+1);
for i = 1 : 2*M-K,
    for j = 1 : K+1,
        S(i, j) = y_hat(K + i - j + 1);
    end
end

%S

% Singular value decomposition
[U, Sigma, V] = svd(S);

% Extract annihilating filter coefficients
a = V(:,K+1);

% Find roots of A(z), the z-transform of the annihilating filter
% Note that multiplying by the correct number of factors of z turns A(z)
% into a polynomial with coefficients a_0, a_1, ... , a_K. Therefore we can
% use matlab's polynomial root function.
u = roots(a);
clear i;
t = mod(real(log(u)*i*tau/(2*pi)), tau);

% Reconstruct spike amplitudes
phi = zeros(N, K);
for i = 1 : N
    for k = 1 : K
        phi(i, k) = sinc_kernel(n(i)-t(k), B, tau);
    end
end
c = pinv(phi)*y;

elapsed_time = toc(ticID);

sigmae = 0;
[t, c] = sort_t_c(t', c');



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

tk
ck
t
c
t_ml
c_ml


% Plot results
%%{
figure;
title('Noiseless samples (z)');
plot(n, z);
xlim([0 tau]);

figure;
title(['Noisy samples (y), SNR = ' num2str(SNR)]);
plot(n, y);
xlim([0 tau]);

figure;
hold on;
stem(tk,ck,'r');
stem(t,c,'b');
stem(t_ml,c_ml,'g');
legend('true', 'cadzow', 'IterML');
xlim([0 tau]);
%%}