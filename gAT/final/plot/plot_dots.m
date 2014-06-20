function h1 = plot_dots(x, y, color, thick)

color = string_to_color(color);

if nargin < 4
   thick = 0;
end

dots = 0;

[style, width, markersize] = get_params(dots, thick);

h1 = plot(x, y, ...
    'Color', color, 'LineStyle', 'none', 'LineWidth', width, 'Marker', 'o', 'MarkerFaceColor', lighten(color), 'MarkerSize', markersize);

end
