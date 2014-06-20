function h1 = plot_line(x, y, color, dots, thick)

% dots - optional param for dotted line (0/1 value)
% thick - optional param for extra thick line (0/1/2 value)

color = string_to_color(color);

if nargin < 4
   dots = 0;
end

if nargin < 5
   thick = 0;
end

[style, width, markersize] = get_params(dots, thick);

h1 = plot(x, y);

set(h1, 'Color', color, 'LineWidth', width, 'LineStyle', style);

end
