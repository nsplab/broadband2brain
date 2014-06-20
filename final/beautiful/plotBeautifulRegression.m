function h = plotBeautifulRegression(X, Y, XRange)

[b, dev, stats] = glmfit(X, Y, 'normal');

b
stats.se

stats = regstats(Y, X, 'linear');
stats.rsquare

hold on;
h = plot(XRange, glmval(b, XRange, 'identity'), '--', 'LineWidth', 3, ...
    'Color', 'black');

% Makes this line not appear in the legend.
set(get(get(h, 'Annotation'), 'LegendInformation'), ...
    'IconDisplayStyle', 'off');
