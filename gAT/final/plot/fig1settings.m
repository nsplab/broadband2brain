function fig = fig1settings

fig = figure;%('units','normalized','position',[.1,.1,.2,.2]);
hold on;
box on;
set(gca, 'XTick', [], 'YTick', []);
xlim([-0.5 8.5]);
ylim([0 10]);

end