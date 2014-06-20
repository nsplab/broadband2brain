% generate a bunch of random signals for use by cramer rao code
% this allows same random signals to be used for each data point (i.e. for
% different values of sigmae)
% t(signal_num, k) = spike time

function [t_mat, c_mat] = random_signals(num, K)

dir = 'C:\Users\alex\Desktop\Sub-Nyquist Signal Processing\figures\cramer_rao\';

t_mat = zeros(num, K);
c_mat = zeros(num, K);

for i = 1 : num

    % generate random spikes
    delta = 2;
    ck = normrnd(10, 4, 1, K);
    %ck = unifrnd(0, 20, 1, K);
    tk = sort(unifrnd(0, 20, 1, K));
    while min(diff(tk)) < delta
        tk = sort(unifrnd(0, 20, 1, K));
    end
    
    t_mat(i, :) = tk;
    c_mat(i, :) = ck;

end

save([dir 'mat/cr_randsigs_' int2str(num)], 't_mat', 'c_mat');

end