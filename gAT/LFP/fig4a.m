% Fig 4a: LFP distribution, all 4 methods
function [] = fig4a()

clear all;
close all;

% Parameters
% Original: elec = 9, T = 0.005
elec = 9;
%T = 0.100;
T = 0.015;
%T = 0.03;
%T = 0.05;
%T = 0.10;

% plot of distributions, all 4 methods
data_set = 1;
methods = {'real', 'AT', 'gAT', 'TD'};
Ts = T*[0 1 1 1];
colors = {[0 0 0], [1 0 0], [0 0 1], [0 1 0]};
%widths = [6 3 3];
%styles = {'--','-','-'};
dots = [1 0 0 0];
thicks = [0 0 0 0];

fig = figure;%('Units','normalize','Position',[1/4 1/4 1/2 1/2]);
hold on;

h = [];

for i = 1 : length(methods)

    load([root_dir() '../LFP/distr/distr_' methods{i} '_' int2str(data_set) '_' int2str(elec) '_' int2str(1000*Ts(i))])

    %subplot(1,2,1), plot(deg,distrlo,'Color',colors{i},'LineWidth',2); hold on; axis square; ylim1 = get(gca,'YLim');
    %p = plot(deg,distrhi,'Color',colors{i},'LineWidth',widths(i),'LineStyle',styles{i});
    p = plot_line(deg, distrhi, colors{i}, dots(i), thicks(i));
    %axis square; ylim2 = get(gca,'YLim');
    h = [h p];
    
    set(gca,'Box','off','XLim',[-180 180]);
    %line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
    xlabel('Phase (degrees)', 'FontSize', 28);
    %xlabel(['Phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz Oscillations (deg)'], 'FontSize', 14);
    ylabel('Relative Risk of Spiking', 'FontSize', 28);% title(['High Amplitude (>' num2str(cutoff,'%5.1f') '\muV)'], 'FontSize', 28);
    %text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nhi))],[' p = ' num2str(phi)],[' pd = ' num2str(pdhi,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
    set(gca, 'FontSize', 24, 'XTick', [-180,-90,0,90,180], 'YTick', [0.6 0.8 1 1.2 1.4]);
    
end

%title(['T = ' int2str(T*1000)], 'FontSize', 28);
ylim([0.6 1.4]);
%ylim([15 35]);
axis manual;
legend(h,'True', 'AT', 'gAT', 'twodelta');
direc = [root_dir() '../paper2plots/final/fig4/'];
saveas(fig, [direc 'fig4a_' int2str(elec) '_' num2str(1000*T)], 'epsc');
saveas(fig, [direc 'fig4a_' int2str(elec) '_'  num2str(1000*T)], 'fig');

fig = figure;%('Units','normalize','Position',[1/4 1/4 1/2 1/2]);
hold on;

h = [];

for i = 1 : length(methods)

    load([root_dir() '../LFP/distr/distr_' methods{i} '_' int2str(data_set) '_' int2str(elec) '_' int2str(1000*Ts(i))]);
    
    % fit cos
    %{
    x = lsqcurvefit(@fun, [0,0,1], deg*pi/180, distrhi');
    p = plot_line(deg, fun(x, deg*pi/180), colors{i}, dots(i), thicks(i));
    %axis square; ylim2 = get(gca,'YLim');
    h = [h p];
    %}
    
    % glmfit instead
    vx = cos(lfp_vec); % not actually velocities but analogous to velocities in tuning curve calculation
    vy = sin(lfp_vec);
    [b, dev, stats] = glmfit([vx vy], spk_vec, 'poisson', 'const', 'on');
    % normalize b_0 = 0
    b(1) = 0;
    p = plot_line(deg, exp(b(1)+b(2)*cos(deg*pi/180)+b(3)*sin(deg*pi/180)), colors{i}, dots(i), thicks(i));
    h = [h p];
    
    set(gca,'Box','off','XLim',[-180 180]);
    %line(get(gca,'XLim'),ones(1,2),'Color','k','LineStyle',':');
    xlabel('Phase (degrees)', 'FontSize', 28);
    %xlabel(['Phase of ' num2str(blfp(1)) '-' num2str(blfp(2)) 'Hz Oscillations (deg)'], 'FontSize', 14);
    ylabel('Relative Risk of Spiking', 'FontSize', 28);% title(['High Amplitude (>' num2str(cutoff,'%5.1f') '\muV)'], 'FontSize', 28);
    %text(min(get(gca,'XLim')),max(get(gca,'YLim')),{[' N = ' num2str(sum(nhi))],[' p = ' num2str(phi)],[' pd = ' num2str(pdhi,'%5.1f') 'deg']},'HorizontalAlignment','left','VerticalAlignment','top');
    set(gca, 'FontSize', 24, 'XTick', [-180,-90,0,90,180], 'YTick', [0.6 0.8 1 1.2 1.4]);
     
    % Plot preferred direction
    %plot_dots(180/pi*x(2), x(1)+x(3), colors{i}, thicks(i));
    b_unit = b(2:3)/norm(b(2:3));
    plot_dots(atan2(b(3),b(2))*180/pi, exp(b(1)+[b(2) b(3)]*b_unit), colors{i}, thicks(i));
    
end

%title(['T = ' int2str(T*1000)], 'FontSize', 28);
ylim([0.6 1.4]);
%ylim([15 35]);
axis manual;
legend(h,'True', 'AT', 'gAT', 'twodelta');
direc = [root_dir() '../paper2plots/final/fig4/'];
saveas(fig, [direc 'figa2_' int2str(elec) '_' num2str(1000*T)], 'epsc');
saveas(fig, [direc 'figa2_' int2str(elec) '_' num2str(1000*T)], 'fig');

end

function y = fun(x, xi)
y = x(1)*cos(xi-x(2))+x(3);
end
