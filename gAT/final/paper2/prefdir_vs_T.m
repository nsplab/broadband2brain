% Plot error in preferred direction vs T for AT and aFRI
% tests whether aFRI gets better tuning curve estimates as T increases
% can also plot tuning curves (using est_delay_b) for particular T (0.015)
% see TESTING area

% TODO
% add AT as a method (DONE)
% AT throw out double-spikes?

%close all;

compute = 0;

if compute

%%{
T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % match uni_params
paramIDs = 0 : 15;

T_vec = T_vec(1:end);
paramIDs = paramIDs(1:end);

method_names = {'AT', 'analog', 'twodelta'};
paramOffset = [10, 10, 10];
data_set = 1;
electrodes = channels_to_use(data_set);

colors = {[0 0 1], [0 1 0], [1 0 1]};

h = [];

data = cell(1, 3);
data{1} = [];  % data{1}(chan)(T_ind) = error
data{2} = [];
data{3} = [];

for i = 1 : length(method_names)
    for j = 1 : length(T_vec)
    %for j = 3 % TESTING: for plotting tuning curves
        %T_vec(j)
        [real_dir_vec dir_vec dir_err dir_err_vec] = est_delay_b(method_names{i}, paramIDs(j)+paramOffset(i), data_set, 1, electrodes);
        data{i}(:, j) = dir_err_vec;
        close all;
    end
end

save([root_dir() '../figures/mat/prefT']);
%%}

end

%%{
%clear all;
load([root_dir() '../figures/mat/prefT']);

dir = [root_dir() '../paper2plots/final/fig5/'];

colors = {[1 0 0], [0 0 1], [0 1 0]};

%close all;
%f = figure;
hold on;

for i = 1 : length(data)
    data{i} = data{i} * 180 / pi;
end

for i = 1 : 3%length(methods)
    ci = zeros(2, length(T_vec));
    mean_vec = zeros(1, length(T_vec));
    for j = 1 : length(T_vec)
        dat = data{i}(:,j);
        ci(:,j) = bootci(1000, @mean, dat);
        mean_vec(j) = mean(dat);
    end
    h1 = plot_errorbar(T_vec*1000, mean_vec, ci(1,:), ci(2,:), colors{i});
    h = [h h1];
end

%leg = legend(h, 'AT', 'gAT', 'twodelta');
%set(leg, 'FontSize', 28, 'Location', 'NorthWest');

xlabel('Sampling Period T (ms)', 'FontSize', 20);
ylabel('Error in Preferred Direction (degrees)', 'FontSize', 20);
ylabel('Preferred Direction Error (degrees)', 'FontSize', 20);
%title(['error in preferred direction vs sampling rate, subject ' int2str(data_set)]);

xlim(1000*[0 0.103]);
ylim([0 80]);

set(gca, 'FontSize', 18, 'XTick', 0:10:100);
%set(gca, 'FontSize', 18, 'XTick', []);
%saveas(f, [dir 'prefdirT'], 'epsc');
%saveas(f, [dir 'prefdirT'], 'fig');
%%}
