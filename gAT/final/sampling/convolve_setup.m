function [] = convolve_setup(T, window, L, h, args)
% Precomputes matrix for convolve_fast, saves it in .mat file
% Saved in .mat file: conv_mat and t_y (sample times relative to t1, start of interval)

filename = [root_dir() 'sampling/conv_mat.mat'];

time = -window + get_dt() : get_dt() : T + window - get_dt();  % Add some buffer to account for numerical bugs
T_s = T / L;
t_y = linspace(T_s, T, L);

conv_mat = h(repmat(t_y', [1 length(time)]) - repmat(time, [L 1]), args) * get_dt();
%conv_mat1 = zeros(L, length(time));
%conv_mat2 = zeros(L, length(time));
%for j = 1 : L
%    for k = 1 : length(time)
%        conv_mat1(j, k) = h(t_y(j) - time(k), args) * get_dt();
%    end
%    conv_mat2(j, :) = h(t_y(j) - time, args) * get_dt();
%end
save(filename, 'conv_mat', 't_y');

end
