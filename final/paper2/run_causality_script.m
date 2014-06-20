% compute granger causality matrix for
% 1) real spikes
% 2) AT paramIDs 10 - 20
% 3) aFRI paramIDs 10 - 20
% results are saved (see test_causality)

% Compute causality matrices

%{

%method = 'real';
%paramID = 0;
%test_causality(method, paramID);

method = 'AT';
for paramID = 10:20
    test_causality(method, paramID);
end

method = 'analog';
for paramID = 10:20
    test_causality(method, paramID);
end

%}


% Plot X

%{

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % from uni_params

method = 'real';
paramID = 0;
[Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
figure;
imagesc(X(1:100,:));
title([method ' ' int2str(paramID)]);

method = 'AT';
for paramID = 10:20
    [Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
    figure;
    imagesc(X(1:100,:));
    title([method ' ' int2str(paramID)]);
end

method = 'analog';
for paramID = 10:20
    [Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
    figure;
    imagesc(X(1:100,:));
    title([method ' ' int2str(paramID)]);
end

%}


% Plot numspikes

%{

close all;

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % from uni_params
electrodes = [3 8 9 10 12 13 14 15 16];  % from test_causality

for j = 1 : 9
    
    numspikes_AT = zeros(1, 10);
    numspikes_aFRI = zeros(1, 10);

    method = 'real';
    paramID = 0;
    [Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
    numspikes_real = nnz(X(:,j));

    method = 'AT';
    for paramID = 10:20
        [Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
        numspikes_AT(paramID-9) = nnz(X(:,j));
    end

    method = 'analog';
    for paramID = 10:20
        [Phi Psi1 Psi2 X] = test_causality(method, paramID, 1);
        numspikes_aFRI(paramID-9) = nnz(X(:,j));
    end
    
    figure;
    hold on;
    plot(T_vec, numspikes_real*ones(size(T_vec)), 'c');
    plot(T_vec, numspikes_AT, 'r');
    plot(T_vec, numspikes_aFRI, 'b');
    legend('real', 'AT', 'aFRI');
    xlabel('sampling period (T)');
    ylabel('number of spikes');
    title(['channel ' int2str(electrodes(j))]);
    ymax = ylim;
    ymax = ymax(2);
    ylim([0 ymax]);
    
end

%}


% Load results and plot matrices

%%{
method = 'real';
paramID = 0;
load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
figure;imagesc(Phi);xlabel('Triggers');ylabel('Targets');title(['Phi ' method ' ' int2str(paramID)]);
figure;imagesc(Psi1);xlabel('Triggers');ylabel('Targets');title(['Psi1 ' method ' ' int2str(paramID)]);
figure;imagesc(Psi2);xlabel('Triggers');ylabel('Targets');title(['Psi2 ' method ' ' int2str(paramID)]);
Psi2

method = 'AT';
for paramID = 10:20
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
    figure;imagesc(Phi);xlabel('Triggers');ylabel('Targets');title(['Phi ' method ' ' int2str(paramID)]);
    figure;imagesc(Psi1);xlabel('Triggers');ylabel('Targets');title(['Psi1 ' method ' ' int2str(paramID)]);
    figure;imagesc(Psi2);xlabel('Triggers');ylabel('Targets');title(['Psi2 ' method ' ' int2str(paramID)]);
end

method = 'analog';
for paramID = 10:20
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
    figure;imagesc(Phi);xlabel('Triggers');ylabel('Targets');title(['Phi ' method ' ' int2str(paramID)]);
    figure;imagesc(Psi1);xlabel('Triggers');ylabel('Targets');title(['Psi1 ' method ' ' int2str(paramID)]);
    figure;imagesc(Psi2);xlabel('Triggers');ylabel('Targets');title(['Psi2 ' method ' ' int2str(paramID)]);
end
%%}


% Plot difference from real as a function of T

% Note: Phi can be changed to Psi1 or Psi2 (but be consistent throughout!)

%{
method = 'real';
paramID = 0;
load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
phi_real = Phi;

param_vec = 10:20;

method = 'AT';
AT_vec = zeros(1, length(param_vec));
for i = 1 : length(param_vec)
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
    AT_vec(i) = norm(Phi-phi_real, 'fro');
end

method = 'analog';
aFRI_vec = zeros(1, length(param_vec));
for i = 1 : length(param_vec)
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
    aFRI_vec(i) = norm(Phi-phi_real, 'fro');
end

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % from uni_params

figure;
hold on;
plot(T_vec, AT_vec, 'r');
plot(T_vec, aFRI_vec, 'b');
legend('AT', 'aFRI');
xlabel('sampling period (T)');
ylabel('Frobenius norm');
title('Phi');
%}


% Submatrix plot (high vs low SNR)

% Note: Phi can be changed to Psi1 or Psi2 (but be consistent throughout!)
% electrodes = [3 8 9 10 12 13 14 15 16];  % from test_causality

%{
high_chan = [1 2 3 8 9];
low_chan = [4 5 6 7];
chan = low_chan;

method = 'real';
paramID = 0;
load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
phi_real = Psi2(chan,chan);

param_vec = 10:20;

method = 'AT';
AT_vec = zeros(1, length(param_vec));
for i = 1 : length(param_vec)
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
    AT_vec(i) = norm(Psi2(chan,chan)-phi_real, 'fro');
end

method = 'analog';
aFRI_vec = zeros(1, length(param_vec));
for i = 1 : length(param_vec)
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
    aFRI_vec(i) = norm(Psi2(chan,chan)-phi_real, 'fro');
end

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % from uni_params

figure;
hold on;
plot(T_vec, AT_vec, 'r');
plot(T_vec, aFRI_vec, 'b');
legend('AT', 'aFRI');
xlabel('sampling period (T)');
ylabel('Frobenius norm');
title('Psi2 (low SNR)');
%}


% Channel-by-channel plot

%{

close all;

T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04 0.05];  % from uni_params
electrodes = [3 8 9 10 12 13 14 15 16];  % from test_causality

for j = 1 : 9

    method = 'real';
    paramID = 0;
    load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
    phi_real = Phi;

    % TESTING: use different "phi_real"
    %method = 'analog';
    %paramID = 10;
    %load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)]);
    %phi_real = Phi;

    param_vec = 10:20;

    method = 'AT';
    AT_vec = zeros(1, length(param_vec));
    for i = 1 : length(param_vec)
        load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
        AT_vec(i) = abs(Phi(j,j)-phi_real(j,j));
    end

    method = 'analog';
    aFRI_vec = zeros(1, length(param_vec));
    for i = 1 : length(param_vec)
        load(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(param_vec(i))]);
        aFRI_vec(i) = abs(Phi(j,j)-phi_real(j,j));
    end

    f = figure;
    hold on;
    plot(T_vec, AT_vec, 'r');
    plot(T_vec, aFRI_vec, 'b');
    legend('AT', 'aFRI');
    xlabel('sampling period (T)');
    ylabel('absolute difference in self-causality');
    title(['channel ' int2str(electrodes(j))]);
    ymax = ylim;
    ymax = ymax(2);
    ylim([0 ymax]);
    print(f, ['C:\Users\alex\Desktop\caus\caus' int2str(j) '.pdf'], '-dpdf', '-r0');

end

%}