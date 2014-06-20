% tests causality between neurons
% adapted from granger causality code
function [Phi Psi1 Psi2 X] = test_causality(method, paramID, Xonly)

% if Xonly = 1 then only return X (return 0 for others)
if nargin < 3
    Xonly = 0;
end

data_set = 1;
% TESTING: eliminated neuron 11 so that X would be the same size as theirs
electrodes = [3 8 9 10 12 13 14 15 16]; %channels_to_use(data_set);
bin_size = 0.001;
num_bins = 100000;
duration = bin_size*num_bins;

if ~strcmp(method, 'real')
    [~, T] = uni_params(method, paramID);
end
    
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start;  % use all available data
start_time = 0.001 * ceil(start_time / 0.001);  % make start time a multiple of 1 ms
end_time = start_time + duration;

time = linspace(start_time, end_time, num_bins);  % length 100000
%median(diff(time))

X = zeros(length(time), length(electrodes));

if strcmp(method, 'real')
    
    for i = 1 : length(electrodes)
        spike_times = struct();
        spike_times(1).data = real_spikes(data_set, electrodes(i), start_time, end_time, 1);
        spike_vec = spike_times_to_vec(time, spike_times')';
        X(:, i) = spike_vec;
    end

else
    
    filename = [root_dir() 'reconstruct/saved_spikes/spk_' method '_' int2str(data_set) '_' int2str(paramID) '.mat'];
    load(filename);  % Loads t, c
    for i = 1 : length(electrodes)
        t_new = t{electrodes(i)}';
        c_new = c{electrodes(i)}';
        ind = find(t_new >= start_time-T & t_new < end_time+T);
        spike_vec = spike_to_vec_an(t_new(ind), c_new(ind), T, bin_size, start_time, num_bins)';
        X(:, i) = spike_vec;
    end

end

%nnz(X == 0)
%nnz(X == 1)
%nnz(X == 2)

if Xonly
    Phi = 0;
    Psi1 = 0;
    Psi2 = 0;
    return
end


% ADAPTED FROM demo_sim.m

%%{

% Load simulated data
%load data_sim_9neuron.mat;     % 9-neuron network
% load data_sim_hidden.mat;      % 5-neuron network with hidden feedback

% Dimension of input data (L: length, N: number of neurons)
[L,N] = size(X);

% To fit GLM models with different history orders
for neuron = 1:N                            % neuron
    for ht = 2:2:10                         % history, when W=2ms
        [bhat{ht,neuron}] = glmwin(X,neuron,ht,200,2);
    end
end

% To select a model order, calculate AIC
for neuron = 1:N
    for ht = 2:2:10
        LLK(ht,neuron) = log_likelihood_win(bhat{ht,neuron},X,ht,neuron,2); % Log-likelihood
        aic(ht,neuron) = -2*LLK(ht,neuron) + 2*(N*ht/2 + 1);                % AIC
    end
end

% To plot AIC
for neuron = 1:N
    figure(neuron);
    plot(aic(2:2:10,neuron));
    title(['neuron ' int2str(neuron)]);
end

% Save results
save(['C:\Users\alex\Desktop\final\paper2\mat\result' method '_' int2str(paramID)],'bhat','aic','LLK');

%%}




%  ADAPTED FROM CausalTest.m

%%{

% Load data
%load data_sim_9neuron.mat;     % 9-neuron network
% load data_sim_hidden.mat;      % 5-neuron network with hidden feedback
%load result_sim.mat;

load(['C:\Users\alex\Desktop\final\paper2\mat\result' method '_' int2str(paramID)]);

% Selected spiking history orders by AIC
%ht = 2*[3 2 3 3 3 2 2 3 3];      % for 9-neuron network
% ht = 2*[5 2 2];                  % for 5-neuron network with hidded feedback

%ht = 2*[5 3 5 1 1 3 4 5 2];% 1];
ht = 2*[2 2 2 2 2 2 2 2 2];

% Dimension of data (L: length, N: number of neurons)
[L,N] = size(X);

% Re-optimizing a model after excluding a trigger neuron's effect and then
% Estimating causality matrices based on the likelihood ratio
for target = 1:N
    LLK0(target) = LLK(ht(target),target);              % Likelihood of full model
    % LLK0(target) = log_likelihood_win(bhat{ht(target),target},X,ht(target),target);
    for trigger = 1:N
        % MLE after excluding trigger neuron
        [bhatc{target,trigger}] = glmcausal(X,target,trigger,ht(target),200,2);
        
        % Log likelihood obtained using a new GLM parameter and data, which exclude trigger
        LLKC(target,trigger) = log_likelihood_causal(bhatc{target,trigger},X,trigger,ht(target),target,2);
        
        % Log likelihood ratio
        LLKR(target,trigger) = LLKC(target,trigger) - LLK0(target);
        
        % Sign (excitation and inhibition) of interaction from trigger to target
        % Averaged influence of the spiking history of trigger on target
        SGN(target,trigger) = sign(sum(bhat{ht(target),target}(ht(target)/2*(trigger-1)+2:ht(target)/2*trigger+1)));
    end
end

% Granger causality matrix, Phi
Phi = -SGN.*LLKR;

% ==== Significance Test ====
% Causal connectivity matrix, Psi, w/o FDR
D = -2*LLKR;                                     % Deviance difference
alpha = 0.05;
for ichannel = 1:N
    temp1(ichannel,:) = D(ichannel,:) > chi2inv(1-alpha,ht(ichannel)/2);
end
Psi1 = SGN.*temp1;

% Causal connectivity matrix, Psi, w/ FDR
fdrv = 0.05; 
temp2 = FDR(D,fdrv,ht);
Psi2 = SGN.*temp2;

% Plot the results
%figure(1);imagesc(Phi);xlabel('Triggers');ylabel('Targets');
%figure(2);imagesc(Psi1);xlabel('Triggers');ylabel('Targets');
%figure(3);imagesc(Psi2);xlabel('Triggers');ylabel('Targets');

% Save results
save(['C:\Users\alex\Desktop\final\paper2\mat\CausalMaps_' method '_' int2str(paramID)],'bhatc','LLK0','LLKC','LLKR','D','SGN','Phi','Psi1','Psi2');

%%}