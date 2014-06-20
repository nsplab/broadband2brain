% Script for producing figure 1
function [] = fig1(data_set, electrode, start_offset, T, num_intervals, thres)

close all;

% Potential data segments:

%electrode = 16;
%start_offset = 1.05;
%T = 0.01;

%electrode = 9;
%start_offset = 0.99;
%T = 0.005;

%electrode = 9;
%start_offset = 0.998;
%T = 0.005;

%electrode = 9;
%start_offset = 0.9975;
%T = 0.003;


% Parameters
%data_set = 1;
%electrode = 9;
%start_offset = 1.0002;
%T = 0.005;
%num_intervals = 3;
%thres = 60;

buffer = 0.001;
duration = T*num_intervals + 2*buffer;

raw_col = [0.5 0.5 0.5];%[21 51 173]/255;%lighten([63 143 210]/255);
hp_col = [21 51 173]/255;
thres_col = [0 191 50]/255;
comp_col = [0.8 0 0];
int_col{1} = [255 124 0]/255;
int_col{2} = [151 2 157]/255;
int_col{3} = [0 255 0]/255;
int_col{4} = [0 0 255]/255;

% Data segment
[training_start, training_end, data_start, data_end] = data_division(data_set);
start_time = data_start + start_offset;
end_time = start_time + duration;

% Get data
buff = 0.1;
[time, dat_raw] = get_data(data_set, electrode, start_time-buff, end_time+buff, 'raw');
time = time - (start_time + buffer);  % start first interval at time 0

% Highpass
[A, B] = butter_bp(300, 6000, 2);  % Change filter order?
args.A_filter = A;
args.B_filter = B;
hp_handle = @filter_generic;  % Use causal version
dat_hp = -hp_handle(dat_raw, args);

% Threshold
%thres = get_thres(data_set, electrode);
%thres = 80;

% Segment intervals
ind = cell(1, num_intervals);
for i = 1 : num_intervals
    ind{i} = find_fast(time, (i-1)*T, i*T, get_dt());
end
ind_all = find_fast(time, 0, num_intervals*T, get_dt());

% Comparator
comp_y = 1;%thres;
comp = comp_y*(dat_hp >= thres);

% Latched comparator
comp_latch = comp;
for i = 1 : num_intervals
    comp_latch(ind{i}) = latch(comp_latch(ind{i}));
end

% 2 successive integrals
x1 = zeros(size(time));
x2 = zeros(size(time));
x3 = zeros(size(time));
x4 = zeros(size(time));
t = zeros(1, num_intervals);
w = zeros(1, num_intervals);
t2 = zeros(2, num_intervals);
w2 = zeros(2, num_intervals);
for i = 1 : num_intervals
    [y xs] = succ_int_2(comp(ind{i}), get_dt(), 4);
    ys{i} = y;
    x1(ind{i}) = xs{1};
    x2(ind{i}) = xs{2};
    x3(ind{i}) = xs{3};
    x4(ind{i}) = xs{4};
    [t(i, :), w(i, :)] = reconstruct_analog(y, T, struct(), []);
    t(i) = t(i) + (i-1)*T;
    [t2(:, i) w2(:, i)] = reconstruct_twodelta(y, T, struct(), []);
    t2(:, i) = t2(:, i) + (i-1)*T;
end

% Plot
f = figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.8]);
xl = [-buffer buffer+num_intervals*T]*1000;

% Raw data (AT)
subplot(5, 6, [2 3]);
hold on;
plot(time*1000, -dat_raw, 'Color', raw_col);
xlim(xl);
ylim([-150 250]);
%xlabel('time (ms)');
%ylabel('voltage (uV)');
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
ylim(yl);

%title(['$\delta y_1 y_2 x_1(t) x_2(t) t_1 t_2$'], 'Interpreter', 'latex', 'FontSize', 14);

% Raw data (aFRI)
%subplot(5, 6, 2);
%hold on;
%plot(time*1000, -dat_raw, 'Color', raw_col);
%xlim(xl);
%ylim([-150 250]);
%%xlabel('time (ms)');
%%ylabel('voltage (uV)');
%box on;
%set(gca, 'xtick', []);
%set(gca, 'ytick', []);
%% T-Intervals
%yl = ylim;
%for i = 0 : num_intervals
%    plot([i*T i*T]*1000, yl, 'k:');
%end
%ylim(yl);

% Raw data (twodelta)
%subplot(5, 6, 3);
%hold on;
%plot(time*1000, -dat_raw, 'Color', raw_col);
%xlim(xl);
%ylim([-150 250]);
%%xlabel('time (ms)');
%%ylabel('voltage (uV)');
%box on;
%set(gca, 'xtick', []);
%set(gca, 'ytick', []);
%% T-Intervals
%yl = ylim;
%for i = 0 : num_intervals
%    plot([i*T i*T]*1000, yl, 'k:');
%end
%ylim(yl);

% Highpassed data (AT)
subplot(5, 6, [8 9]);
hold on;
plot(time*1000, dat_hp, 'Color', hp_col);
xlim(xl);
ylim([-170 220]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
%edge_indices = ind_all([comp(ind_all(1:(end-1))) ~= comp(ind_all(2:end))]);
%edges = time(edge_indices) + median(diff(time)) / 2;
up_edge_indices = ind_all([comp(ind_all(1:(end-1))) == 0 & comp(ind_all(2:end)) == 1]);
down_edge_indices = ind_all([comp(ind_all(1:(end-1))) == 1 & comp(ind_all(2:end)) == 0]);
edge_indices = [up_edge_indices (down_edge_indices + 1)];
%edges = time(edge_indices);
edges = time(edge_indices) + median(diff(time));
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, [yl(1) dat_hp(edge_indices(i))], 'k:');
end
ylim(yl);
% Threshold
plot(xl, [thres thres], 'Color', thres_col, 'LineStyle', '--', 'LineWidth', 2);


% Highpassed data (aFRI)
%subplot(5, 6, 5);
%hold on;
%plot(time*1000, dat_hp, 'Color', hp_col);
%xlim(xl);
%ylim([-170 220]);
%box on;
%set(gca, 'xtick', []);
%set(gca, 'ytick', []);
%% T-Intervals
%yl = ylim;
%for i = 0 : num_intervals
%    plot([i*T i*T]*1000, yl, 'k:');
%end
%ylim(yl);
%% Threshold
%plot(xl, [thres thres], 'Color', thres_col, 'LineStyle', '--', 'LineWidth', 2);

% Highpassed data (twodelta)
%subplot(5, 6, 6);
%hold on;
%plot(time*1000, dat_hp, 'Color', hp_col);
%xlim(xl);
%ylim([-170 220]);
%box on;
%set(gca, 'xtick', []);
%set(gca, 'ytick', []);
%% T-Intervals
%yl = ylim;
%for i = 0 : num_intervals
%    plot([i*T i*T]*1000, yl, 'k:');
%end
%ylim(yl);
%% Threshold
%plot(xl, [thres thres], 'Color', thres_col, 'LineStyle', '--', 'LineWidth', 2);

% Latched comparator (AT)
subplot(5, 6, [13 14]);
hold on;
for i = 1 : num_intervals
    plot_box(time(ind{i})*1000, comp_latch(ind{i}), comp_col);
end
xlim(xl);
ylim([-0.1 1.1]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', [0 1]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, yl, 'k:');
end
ylim(yl);

% Comparator (aFRI)
subplot(5, 6, [15 16]);
hold on;
plot_box(time(ind_all)*1000, comp(ind_all), comp_col);
xlim(xl);
ylim([-0.1 1.1]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', [0 1]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, yl, 'k:');
end
ylim(yl);

% Comparator (twodelta)
%subplot(5, 6, 9);
%hold on;
%plot_box(time(ind_all)*1000, comp(ind_all), comp_col);
%xlim(xl);
%ylim([-0.1 1.1]);
%box on;
%set(gca, 'xtick', []);
%set(gca, 'ytick', [0 1]);
%% T-Intervals
%yl = ylim;
%for i = 0 : num_intervals
%    plot([i*T i*T]*1000, yl, 'k:');
%end
%ylim(yl);

% 0/1 (AT)
subplot(5, 6, [19 20]);
hold on;
for i = 1 : num_intervals
    plot_box(time(ind{i})*1000, comp_latch(ind{i}), comp_col);
    i*T*1000
    comp_latch(ind{i}(end))
    plot(i*T*1000, comp_latch(ind{i}(end)), 'Marker', 'o', 'MarkerEdgeColor', comp_col, 'MarkerFaceColor', lighten(comp_col), 'LineWidth', 2);
end
xlim(xl);
ylim([-0.1 1.1]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', [0 1]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, yl, 'k:');
end
ylim(yl);

%factor = [1 300 300^2 300^3];
ratio = 365;
factor = [1 ratio ratio^2 ratio^3];

% Integrals (aFRI)
subplot(5, 6, [21 22]);
hold on;
for i = 1 : num_intervals
    plot_int(time(ind{i})*1000, 10^4*x1(ind{i})*factor(1), 10^4*x2(ind{i})*factor(2), int_col{1}, int_col{2});
    plot(i*T*1000, 10^4*ys{i}(1)*factor(1), 'Marker', 'o', 'MarkerEdgeColor', int_col{1}, 'MarkerFaceColor', lighten(int_col{1}), 'LineWidth', 2);
    plot(i*T*1000, 10^4*ys{i}(2)*factor(2), 'Marker', 'o', 'MarkerEdgeColor', int_col{2}, 'MarkerFaceColor', lighten(int_col{2}), 'LineWidth', 2);
end
xlim(xl);
ylim([-0.5 4.3]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', [0]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, yl, 'k:');
end
ylim(yl);

% Integrals (twodelta)
subplot(5, 6, [23 24]);
hold on;
for i = 1 : num_intervals
    plot_int(time(ind{i})*1000, 10^4*x1(ind{i})*factor(1), 10^4*x2(ind{i})*factor(2), int_col{1}, int_col{2}); % Actual integral
    plot_int(time(ind{i})*1000, 10^4*x3(ind{i})*factor(3), 10^4*x4(ind{i})*factor(4), int_col{3}, int_col{4}); % Actual integral
    %plot_int(time(ind{i})*1000, 10^4*x1(ind{i}), 10^4*x3(ind{i})*300, 10^4*x4(ind{i})*300, 10^4*x2(ind{i})*300, int_col{1}, int_col{2}, int_col{3}, int_col{4}); % Actual integral
    plot(i*T*1000, 10^4*ys{i}(1)*factor(1), 'Marker', 'o', 'MarkerEdgeColor', int_col{1}, 'MarkerFaceColor', lighten(int_col{1}), 'LineWidth', 2); % Filled circle for first
    plot(i*T*1000, 10^4*ys{i}(2)*factor(2), 'Marker', 'o', 'MarkerEdgeColor', int_col{2}, 'MarkerFaceColor', lighten(int_col{2}), 'LineWidth', 2); % Filled circle for first
    plot(i*T*1000, 10^4*ys{i}(3)*factor(3), 'Marker', 'o', 'MarkerEdgeColor', int_col{3}, 'MarkerFaceColor', lighten(int_col{3}), 'LineWidth', 2); % Filled circle for first
    plot(i*T*1000, 10^4*ys{i}(4)*factor(4), 'Marker', 'o', 'MarkerEdgeColor', int_col{4}, 'MarkerFaceColor', lighten(int_col{4}), 'LineWidth', 2); % Filled circle for first
end
for i = 1:size(edges, 2)
    plot([edges(i) edges(i)]*1000, yl, 'k:');
end
xlim(xl);
ylim([-0.5 4.3]);
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', [0]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
ylim(yl);

% Reconstructed (AT)
subplot(5, 6, [25 26]);
hold on;
plot(time*1000, dat_hp, 'Color', lighten(hp_col));
for i = 1 : num_intervals
    if (comp_latch(ind{i}(end)) == 1)
        plot([(i-1)*T+T/2 (i-1)*T+T/2]*1000, ylim, 'Color', 'k', 'LineWidth', 2);
    end
end
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
xlim(xl);
ylim([-170 220]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
ylim(yl);

% Reconstructed (aFRI)
subplot(5, 6, [27 28]);
hold on;
plot(time*1000, dat_hp, 'Color', lighten(hp_col));
for i = 1 : num_intervals
    if w(i) > 0
        plot([t(i) t(i)]*1000, ylim, 'Color', 'k', 'LineWidth', 2);
        %plot_box([t(i)-w(i)/2 t(i)-w(i)/2 t(i)+w(i)/2]*1000, [0 comp_y 0], 'k');
    end
end
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
xlim(xl);
ylim([-170 220]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
ylim(yl);

% Reconstructed (twodelta)
subplot(5, 6, [29 30]);
hold on;
plot(time*1000, dat_hp, 'Color', lighten(hp_col));
for i = 1 : num_intervals
    if w2(1, i) > 0
        plot([t2(1, i) t2(1, i)]*1000, ylim, 'Color', 'k', 'LineWidth', 2);
    end
    if w2(2, i) > 0
        plot([t2(2, i) t2(2, i)]*1000, ylim, 'Color', 'k', 'LineWidth', 2);
    end
end
box on;
set(gca, 'xtick', []);
set(gca, 'ytick', []);
xlim(xl);
ylim([-170 220]);
% T-Intervals
yl = ylim;
for i = 0 : num_intervals
    plot([i*T i*T]*1000, yl, 'k:');
end
ylim(yl);

saveas(f, [root_dir() '../paper2plots/final/fig1'], 'epsc');

end

function output = latch(data)
val = 0;
for i = 1 : length(data)
    data(i) = max(data(i), val);
    val = max(val, data(i));
end
output = data;
end

% Plot square wave with true vertical edges
function h = plot_box(time, dat, col)
last_x = time(1);
last_y = dat(1);
for i = 2 : length(time)
    if dat(i) ~= last_y || i == length(time)
        plot([last_x, time(i)], [last_y last_y], 'Color', col);
        h = plot([time(i) time(i)], [last_y dat(i)], 'Color', col);
        last_x = time(i);
        last_y = dat(i);
    end
end
end

% Plot integrals such that you can see both if they're both 0 (dotted)
% Assumes a leading segment of 0 followed by a segment of nonzero
function [h1 h2] = plot_int(time, x1, x2, c1, c2)
zero_ind = find(x1 == 0);
plot(time(zero_ind), x1(zero_ind), 'Color', c1, 'LineWidth', 2);
plot(time(zero_ind), x2(zero_ind), 'Color', c2, 'LineWidth', 2, 'LineStyle', '--');
nz_ind = zero_ind(end) : length(time);
h1 = plot(time(nz_ind), x1(nz_ind), 'Color', c1, 'LineWidth', 2);
h2 = plot(time(nz_ind), x2(nz_ind), 'Color', c2, 'LineWidth', 2);
end

% Plot integrals such that you can see all four if they're both 0 (dotted)
% Assumes a leading segment of 0 followed by a segment of nonzero
function [h1 h2 h3 h4] = plot_int4(time, x1, x2, c1, c2)
zero_ind = find(x1 == 0);
plot(time(zero_ind), x1(zero_ind), 'Color', c1, 'LineWidth', 2);
plot(time(zero_ind), x2(zero_ind), 'Color', c2, 'LineWidth', 2, 'LineStyle', '--');
nz_ind = zero_ind(end) : length(time);
h1 = plot(time(nz_ind), x1(nz_ind), 'Color', c1, 'LineWidth', 2);
h2 = plot(time(nz_ind), x2(nz_ind), 'Color', c2, 'LineWidth', 2, 'LineStyle', '--');
end
