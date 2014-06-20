function [b_0, b_1, b_2] = est_b(data_set, id, load_flag, spike_vector, vx, vy)
% Estimates parameters b_0, b_1, b_2 using glmfit
% Or just loads them from a .mat file if load_flag = 1
% Assumes 1 neuron per channel

% Inputs:
% data_set
% id - identifier string, allows for different sets of parameters etc. (name of object saved to file)
% load_flag - 0 or 1
% spike_vector(neuron_num, :) = 0's and 1's for presence of spikes

% Outputs:
% each is a row vector with one entry per neuron

dt = get_dt_info();

filename = [root_dir() 'pp_decode/b.mat'];
load(filename);

if load_flag == 1,
    b = eval([id '(data_set)']);
    b_0 = b.b_0;
    b_1 = b.b_1;
    b_2 = b.b_2;
else
    b_0 = zeros(1, C);
    b_1 = zeros(1, C);
    b_2 = zeros(1, C);
    elecs = channels_to_use(data_set);
    assert length(elecs) == size(spike_vector, 1);
    for i = 1 : length(elecs),
        b = glmfit([vx, vy], spike_vector(i, :), 'poisson', 'const', 'on');
        b_0(i) = b(1) - log(dt);  % Correct for the fact that glmfit's lambda is spikes per timestep (not per second)
        b_1(i) = b(2);
        b_2(i) = b(3);
    end
    obj.b_0 = b_0;
    obj.b_1 = b_1;
    obj.b_2 = b_2;
    eval([id '(data_set) = obj']);
    save(filename, id, '-append');
end
