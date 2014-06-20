% Implementation of Blu "Sparse Sampling of Signal Innovations" with sinc
% kernel

function [t c elapsed_time] = run_cadzow_sinc(tau, K, N, SNR, tk, ck)

% Parameters
%{
tau = 0.5;  % [0, tau] is time window containing spikes (20)
K = 7;  % number of spikes (5)
N = 71;  % number of samples (50)
%sigmae = 0;  % (2)
SNR = 5;
%}

T_s = tau/N;
Btau = N;  % odd integer (best when N or N-1)
if mod(N,2) == 0
    Btau = N-1;
end
B = Btau/tau;
M = floor(B*tau/2);

n = linspace(T_s, tau, N);  % points at which samples are taken

% Take samples (sinc kernel)
z = zeros(length(n),1);
for k = 1:length(tk)
    z = z + ck(k)*sinc_kernel(n-tk(k), B, tau)';
end

% Calc sigmae to achieve desired SNR
sigmae = sqrt(sum(z.^2)/(N*10^(SNR/10)));
%SNR_cad = 10*log(sum(z.^2)/(N*sigmae^2))/log(10)

% Add noise
e = normrnd(0, sigmae, length(z), 1);
y = z + e;

%y

%%% Beginning of reconstruction algorithm

% Compute DFT y_hat (can apparently be made faster using FFT)
y_hat = zeros(1, Btau);
clear i;
for m = -M : M
    for nn = 1 : N
        y_hat(m+M+1) = y_hat(m+M+1) + y(nn)*exp(-i*2*pi*nn*m/N);
    end
end

%y_hat


% Start timing
ticID = tic;  % TIME FFT (plus the rest)
y_fft = fft(y, Btau)';

% Cadzow denoising step
%dn_ID = tic;
y_hat = cadzow_denoise(y_hat, K);
%denoise_time = toc(dn_ID)

%y_hat

% Form matrix S
S = zeros(2*M-K, K+1);
for p = 1 : 2*M-K,
    for j = 1 : K+1,
        S(p, j) = y_hat(K + p - j + 1);
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
t = zeros(K, 1);
t(1 : length(u)) = mod(real(log(u)*i*tau/(2*pi)), tau);
if length(u) < K
    disp('WARNING: Cadzow: not enough roots found');
end

% Reconstruct spike amplitudes
phi = zeros(N, K);
for p = 1 : N
    for k = 1 : K
        phi(p, k) = sinc_kernel(n(p)-t(k), B, tau);
    end
end
c = pinv(phi)*y;

elapsed_time = toc(ticID);

sigmae = 0;
[t, c] = sort_t_c(t', c');

% Plot results
%{
figure;
title('Noiseless samples (z)');
plot(z);

figure;
title(['Noisy samples (y), SNR = ' num2str(SNR)]);
plot(y);

figure;
hold on;
stem(tk,ck,'r');
stem(t,c,'b');
legend('original', 'estimated');
%}