function h1 = plot_stem(x, y, color, dots, thick, base)

color = string_to_color(color);

if nargin < 4
   dots = 0;
end

if nargin < 5
   thick = 0;
end

if nargin < 6
    base = zeros(size(x));
end

[style, width, markersize] = get_params(dots, thick);

hold on;
for i = 1 : length(x)
    h1 = plot([x(i) x(i)], [base(i) y(i)]);
    set(h1, 'LineStyle', style, 'Color', color, 'LineWidth', ceil(width/2), 'Marker', 'none');
end

h2 = plot(x, y);
set(h2, 'LineStyle', 'none', 'Color', color, 'LineWidth', width, 'Marker', 'o', 'MarkerFaceColor', lighten(color), 'MarkerSize', markersize);

end
