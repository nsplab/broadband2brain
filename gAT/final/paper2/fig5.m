% Fig 5: example tuning curves, all 3 methods
% Loads data saved by est_deltay_b, which is in turn run by prefdir_vs_T
% prefdir_vs_T makes fig 5a

function [] = fig5(data_set, elec, paramID, plot_raw, plot_fit)

dir = [root_dir() 'paper2/mat/tuning'];

%data_set = 1;
%elec = 11;
%paramID = 25;

method_name = 'AT';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);

y_upper = max([freq_real/get_dt_info() freq_real_b]);
y_lower = min([freq_real/get_dt_info() freq_real_b]);

methods = {'AT', 'analog', 'twodelta'};
for i = 1:size(methods, 2)
    method_name = methods{i};
    load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
    y_upper = max([y_upper freq/get_dt_info() freq_b]);
    y_lower = min([y_lower freq/get_dt_info() freq_b]);
end

if (nargin <= 3 || plot_raw)

method_name = 'AT';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);


% RAW CURVES
%fig = figure;
hold on;

% real spikes raw
p1 = plot_line(180/pi*theta_real, freq_real / get_dt_info(), 'k', 1);

% AT raw
p2 = plot_line(180/pi*theta, freq / get_dt_info(), 'r');

% gAT raw
method_name = 'analog';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
p3 = plot_line(180/pi*theta, freq / get_dt_info(), 'b');
temp = freq;

% TD raw
method_name = 'twodelta';
[dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
p4 = plot_line(180/pi*theta, freq / get_dt_info(), 'g');

xlabel('Angle (degrees)', 'FontSize', 20);
ylabel('Av. Firing Rate (spikes/sec)', 'FontSize', 20);
ylabel('Firing Rate (spikes/sec)', 'FontSize', 20);
xlim([-180 180]);
ylim([0 y_upper]);
%leg = legend([p1 p2 p3 p4], 'True', 'AT', 'gAT', 'twodelta');
%set(leg, 'FontSize', 18, 'Location', 'NorthWest');
set(gca, 'FontSize', 20, 'XTick', [-180,-90,0,90,180]);

%save_dir = [root_dir() '../paper2plots/final/fig5/'];
%saveas(fig, [save_dir 'tuning_raw_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)], 'epsc');
%saveas(fig, [save_dir 'tuning_raw_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)], 'fig');

end




if (nargin <= 4 || plot_fit)

% FIT CURVES B
%fig = figure;
hold on;

% real spikes b
p1 = plot_line(180/pi*theta_real, freq_real_b, [0 0 0], 1);
[val ind] = max(freq_real_b);
plot_dots(180/pi*pref_dir_real, freq_real_b(ind), [0 0 0]);

% AT b
method_name = 'AT';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
p2 = plot_line(180/pi*theta, freq_b, [1 0 0]);
[val ind] = max(freq_b);
plot_dots(180/pi*pref_dir, freq_b(ind), [1 0 0]);

% gAT b
method_name = 'analog';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
p3 = plot_line(180/pi*theta, freq_b, [0 0 1]);
[val ind] = max(freq_b);
plot_dots(180/pi*pref_dir, freq_b(ind), [0 0 1]);

% twodelta c
method_name = 'twodelta';
load([dir '/tuning_' method_name '_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)]);
p4 = plot_line(180/pi*theta, freq_b, [0 1 0]);
[val ind] = max(freq_b);
plot_dots(180/pi*pref_dir, freq_b(ind), [0 1 0]);

xlabel('Angle (degrees)', 'FontSize', 20);
ylabel('Av. Firing Rate (spikes/sec)', 'FontSize', 20);
ylabel('Firing Rate (spikes/sec)', 'FontSize', 20);
xlim([-180 180]);
ylim([0 y_upper]);
%leg = legend([p1 p2 p3 p4], 'True', 'AT', 'gAT', 'twodelta');
%set(leg, 'FontSize', 18, 'Location', 'NorthWest');
set(gca, 'FontSize', 20, 'XTick', [-180,-90,0,90,180]);

%save_dir = [root_dir() '../paper2plots/final/fig5/'];
%saveas(fig, [save_dir 'tuning_b_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)], 'epsc');
%saveas(fig, [save_dir 'tuning_b_' int2str(data_set) '_' int2str(elec) '_' int2str(paramID)], 'fig');

%{
plot_line(theta_real, freq_real / get_dt_info(), 'k');  % Normalize lambda correctly            
p1 = plot_line(180/pi*theta_real, freq_real_b, [0 0 0]);
[val ind] = max(freq_real_b);
plot_dots(180/pi*pref_dir_real, freq_real_b(ind), [0 0 0]);
p2 = plot_line(180/pi*theta, freq_b, [0 0 1]);
[val ind] = max(freq_b);
plot_dots(180/pi*pref_dir, freq_b(ind), [0 0 1]);
xlabel('Angle (degrees)', 'FontSize', 18);
ylabel('Av. Firing Rate (spikes/sec)', 'FontSize', 18);
axis tight;
leg = legend([p1 p2 p3], 'True', 'AT', 'gAT');
set(leg, 'FontSize', 18, 'Location', 'NorthWest');
set(gca, 'FontSize', 14);
%}

end

end
