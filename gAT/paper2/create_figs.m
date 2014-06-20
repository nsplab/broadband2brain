clear;
close all;




h = figure;


subplot(3, 3, [1 2 3]);
prefdir_vs_T;

top = 35;

subplot(3, 3, 4);
fig5(1, 3, 15, 0, 1);
ylim([0 top]);

subplot(3, 3, 5);
fig5(1, 3, 20, 0, 1);
ylim([0 top]);

subplot(3, 3, 6);
fig5(1, 3, 25, 0, 1);
ylim([0 top]);

subplot(3, 3, 7);
fig5(1, 15, 15, 0, 1);
ylim([0 top]);

subplot(3, 3, 8);
fig5(1, 15, 20, 0, 1);
ylim([0 top]);

subplot(3, 3, 9);
fig5(1, 15, 25, 0, 1);
ylim([0 top]);

%p = [ 6.5000
%     13.5000
%      8.0000
%     82.0000
%      4.4000
%      2.4000
%      2.5600
%      1.8000
%      0.6800];
%
%e = [ 96
%     100
%     100
%      32
%     128
%     128
%      64
%      16
%      16];
%
%h = figure;
%hold on;
%scatter(e, p, 'Black', 'Fill');
%
%best = min(p ./ e);
%xl = xlim;
%plot(xl, best * xl, 'Black', 'LineWidth', 2);
%
%thres = (6.3*3.2) / 24 / 365.25 / 10 * 1000;
%plot(xl, [thres thres], 'LineWidth', 2);
%
%xlabel('Simultaneously Recorded Neurons');
%ylabel('Power Usage (mW)');
%
%ylim([0 8]);
%
%saveas(h, [root_dir() '../paper2/power'], 'epsc');
%saveas(h, [root_dir() '../paper2/power'], 'fig');




%h1 = openfig([root_dir() '../paper2plots/final/analyze/n_spike_interval.fig']);
%ax1 = gca; % get handle to axes of figure
%fig1 = get(ax1, 'children'); % get handle to all the children in the figure
%
%h2 = openfig([root_dir() '../paper2plots/final/fig2/frac_mean_all.fig']);
%ax2 = gca; % get handle to axes of figure
%fig2 = get(ax2, 'children'); % get handle to all the children in the figure
%
%h3 = openfig([root_dir() '../paper2plots/final/fig2/t_error_mean_all.fig']);
%ax3 = gca; % get handle to axes of figure
%fig3 = get(ax3, 'children'); % get handle to all the children in the figure
%
%h = figure;
%
%fs = 20;
%
%s1 = subplot(6, 1, [1 2]);
%xlabel('Sampling Period (ms)', 'FontSize', fs);
%ylabel('Fraction of Intervals', 'FontSize', fs);
%
%s2 = subplot(6, 1, [5 6]);
%xlabel('Sampling Period (ms)', 'FontSize', fs);
%ylabel('Valid Fraction', 'FontSize', fs);
%
%s3 = subplot(6, 1, [3 4]);
%xlabel('Sampling Period (ms)', 'FontSize', fs);
%ylabel('Average Spike Time Error (ms)', 'FontSize', fs);
%ylabel('Spike Time Error (ms)', 'FontSize', fs);
%
%
%fig1 = copyobj(fig1, s1); % copy children to new parent aces i.e. the subplot axes
%s1 = subplot(6, 1, [1 2]);
%%legend(fig1(end:-2:2), '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Location', 'EastOutside');
%legend(fig1(end:-2:2), '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'Location', 'Best');
%xlim([0 106]);
%set(gca, 'xtick', 0:10:100);
%set(gca, 'FontSize', fs);
%
%fig2 = copyobj(fig2, s2); % copy children to new parent aces i.e. the subplot axes
%s2 = subplot(6, 1, [5 6]);
%legend(fig2([66 33]), 'AT / gAT-1', 'gAT-2', 'Location', 'Best');
%xlim([0 106]);
%set(gca, 'xtick', []);
%set(gca, 'xtick', 0:10:100);
%set(gca, 'FontSize', fs);
%
%fig3 = copyobj(fig3, s3); % copy children to new parent aces i.e. the subplot axes
%s3 = subplot(6, 1, [3 4]);
%legend(fig3([100 99 66 33]), 'AT Theoretical', 'AT', 'gAT-1', 'gAT-2', 'Location', 'Best');
%xlim([0 106]);
%set(gca, 'xtick', []);
%set(gca, 'xtick', 0:10:100);
%set(gca, 'FontSize', fs);

%saveas(h, [root_dir() '../paper2/basic_comparison'], 'epsc');
%saveas(h, [root_dir() '../paper2/basic_comparison'], 'fig');






%clear;
%
%h(1) = openfig([root_dir() '../paper2plots/final/fig3/ROC_median_15.fig']);
%ax(1) = gca; % get handle to axes of figure
%fig{1} = get(ax(1), 'children'); % get handle to all the children in the figure
%
%h(2) = openfig([root_dir() '../paper2plots/final/fig3/ROC_median_50.fig']);
%ax(2) = gca; % get handle to axes of figure
%fig{2} = get(ax(2), 'children'); % get handle to all the children in the figure
%
%h(3) = openfig([root_dir() '../paper2plots/final/fig3/ROC_median_100.fig']);
%ax(3) = gca; % get handle to axes of figure
%fig{3} = get(ax(3), 'children'); % get handle to all the children in the figure
%
%h(4) = openfig([root_dir() '../paper2plots/final/fig4/figa2_9_15.fig']);
%ax(4) = gca; % get handle to axes of figure
%fig{4} = get(ax(4), 'children'); % get handle to all the children in the figure
%
%h(5) = openfig([root_dir() '../paper2plots/final/fig4/figa2_9_50.fig']);
%ax(5) = gca; % get handle to axes of figure
%fig{5} = get(ax(5), 'children'); % get handle to all the children in the figure
%
%h(6) = openfig([root_dir() '../paper2plots/final/fig4/figa2_9_100.fig']);
%ax(6) = gca; % get handle to axes of figure
%fig{6} = get(ax(6), 'children'); % get handle to all the children in the figure
%
%h = figure;
%for i = 1:3
%  s(i) = subplot(4, 6, 2*(i-1) + [1 2 7 8]);
%  fig{i} = copyobj(fig{i}, s(i));
%end
%for i = 4:6
%  s(i) = subplot(4, 6, 2*(i-4) + [13 14 19 20]);
%  fig{i} = copyobj(fig{i}, s(i));
%end
%
%%saveas(h, [root_dir() '../paper2/roc_LFP'], 'epsc');
%fs = 20;
%
%s(1) = subplot(4, 6, [1 2 7 8]);
%size(fig{1})
%legend(fig{1}(end-1:-2:end-5), 'AT', 'gAT-1', 'gAT-2', 'Location', 'Best');
%s(4) = subplot(4, 6, [1 2 7 8]);
%legend(fig{4}(end:-2:1), 'True', 'AT', 'gAT-1', 'gAT-2', 'Location', 'Best');
%
%for i = 1:3
%    s(i) = subplot(4, 6, 2*(i-1) + [1 2 7 8]);
%    axis([0 1 0 1]);
%    xlabel('False Negatives per True Spike', 'FontSize', fs);
%    ylabel('False Positives per True Spike', 'FontSize', fs);
%    set(gca, 'FontSize', fs);
%end
%
%for i = 4:6
%    s(i) = subplot(4, 6, 2*(i-4) + [13 14 19 20]);
%    xlim([-180 180]);
%    xlabel('Phase (degrees)', 'FontSize', fs);
%    ylabel('Relative Risk of Spiking', 'FontSize', fs);
%    set(gca, 'XTick', [-180 -90 0 90 180]);
%    set(gca, 'FontSize', fs);
%end

%saveas(h, [root_dir() '../paper2/roc_LFP'], 'fig');


%fs = 24;
%
%
%
%
%h(3) = openfig([root_dir() '../paper2plots/final/fig4/fig4b.fig']);
%ax(3) = gca; % get handle to axes of figure
%fig{3} = get(ax(3), 'children'); % get handle to all the children in the figure
%
%h(4) = openfig([root_dir() '../paper2plots/final/fig4/fig4c.fig']);
%ax(4) = gca; % get handle to axes of figure
%fig{4} = get(ax(4), 'children'); % get handle to all the children in the figure
%
%h = figure;
%
%subplot(2, 2, 1);
%hold on;
%p = -180:0.01:180;
%plot_line(p, 1+0.5*cos(p/180*pi), [0 0 0], 1, 0);
%plot_dots(0, 1.5, [0 0 0], 0);
%plot_line(p, 1+0.35*cos(p/180*pi+0.4*pi), [0 0.5 0.5], 0, 0);
%plot_dots(-180*0.4, 1.35, [0 0.5 0.5], 0);
%xlim([-180 180]);
%xlabel('Phase (degrees)', 'FontSize', fs);
%ylabel('Relative Risk of Spiking', 'FontSize', fs);
%set(gca, 'XTick', [-180 -90 0 90 180]);
%set(gca, 'FontSize', fs);
%
%subplot(2, 2, 2);
%hold on;
%p = -180:0.01:180;
%p1 = plot_line(p, 1+0.5*cos(p/180*pi), [0 0 0], 1, 0);
%plot_dots(0, 1.5, [0 0 0], 0);
%p2 = plot_line(p, 1+0.35*cos(p/180*pi+0.4*pi), [0 0.5 0.5], 0, 0);
%plot_dots(-180*0.4, 1.35, [0 0.5 0.5], 0);
%xlim([-180 180]);
%xlabel('Phase (degrees)', 'FontSize', fs);
%ylabel('Relative Risk of Spiking', 'FontSize', fs);
%set(gca, 'XTick', [-180 -90 0 90 180]);
%legend([p1 p2], 'True', 'Estimated');
%set(gca, 'FontSize', fs);
%
%
%for i = 3:4
%    s(i) = subplot(2, 2, i);
%    fig{i} = copyobj(fig{i}, s(i));
%    xlim(1000*[0 0.103]);
%    xlabel('Sampling Period (ms)', 'FontSize', fs);
%    set(gca, 'FontSize', fs);
%end
%
%subplot(2, 2, 3);
%ylabel('Preferred Phase Error (degrees)', 'FontSize', fs);
%
%subplot(2, 2, 4);
%ylabel('Average Distribution Error', 'FontSize', fs);
%legend(fig{3}([99 66 33]), 'AT', 'gAT-1', 'gAT-2');

%saveas(h, [root_dir() '../paper2/LFP'], 'epsc');
%saveas(h, [root_dir() '../paper2/LFP'], 'fig');
