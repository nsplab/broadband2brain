function [y, sample_times] = convolve_fast(x, t1)
% Takes samples by convolution using matrix precomputed by convolve_setup
% x should be a (column) data vector for time [t1 - window, t1 + T + window)

% NOTE: need to run convolve_setup first

filename = [root_dir() 'sampling/conv_mat.mat'];
load(filename);  % Loads conv_mat and t_y

len_x = size(conv_mat, 2);  % Should be length of x
offset = round((length(x) - len_x) / 2);

if length(x) < len_x

    disp('WARNING: length(x) < len_x, convolve_fast.m');
    length_x = length(x)
    len_x
    if size(x, 1) < size(x, 2)
        x = x';
    end
    y = conv_mat(:, 1 : length(x)) * x;

else

    y = conv_mat * x(offset + 1 : offset + len_x)';

end

sample_times = t_y + t1;

end
