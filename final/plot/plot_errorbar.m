function h1 = plot_errorbar(x, y, l, h, color, barsoff, width)

if nargin < 6
    barsoff = 0;
end

show_thres = -1;  % if errorbars are smaller than this, don't show them
% samples: 0.006
% iter: 0.0025

color = string_to_color(color);

plot(x, y, ...
    'Color', lighten(color), 'LineWidth', 3);

for i = 1 : length(x)
    
    if ~barsoff && max(y(i)-l(i), h(i)-y(i)) > show_thres
    
        h1 = errorbar(x(i), y(i), ...
            y(i)-l(i), h(i)-y(i), ...
            'Color', color, 'LineWidth', 3, 'Marker', 'o', 'MarkerFaceColor', lighten(color), 'MarkerSize', 8);

        bar_width = 30;  % normal:30, zoom:60
        if nargin >= 7
            bar_width = width;
        end
        errorbar_tick(h1, bar_width);  % width of bars
        
    end
        
    h1 = plot(x(i), y(i), ...
        'Color', color, 'LineWidth', 3, 'Marker', 'o', 'MarkerFaceColor', lighten(color), 'MarkerSize', 8);

end

end
