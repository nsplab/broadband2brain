function [] = plotBeautifulData(dataMeans, dataLowerBounds, ...
    dataUpperBounds, parSetColor, location, startTrial, fill, legendEntry)

if(nargin < 6)
    startTrial = 1;
    fill = 'filled';
    legendEntry = 'line';
elseif(nargin < 7)
    fill = 'filled';
    legendEntry = 'line';
elseif(nargin < 8)
    legendEntry = 'line';
end

hold on;
h1 = errorbar(startTrial:(length(dataMeans) + startTrial - 1), dataMeans, ...
    dataMeans - dataLowerBounds, dataUpperBounds - dataMeans, ...
    'Color', parSetColor, 'LineWidth', 3);
h2 = [];

if(strcmp(fill, 'filled'))
    h2 = plot(startTrial:(length(dataMeans) + startTrial - 1), dataMeans, ...
        'Color', parSetColor, 'LineStyle', 'none', 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerFaceColor', lightenColor(parSetColor), ...
        'MarkerSize', 8);
elseif(strcmp(fill, 'empty'))
    h2 = plot(startTrial:(length(dataMeans) + startTrial - 1), dataMeans, ...
        'Color', parSetColor, 'LineStyle', 'none', 'LineWidth', 2, ...
        'Marker', 'o', 'MarkerFaceColor', [1 1 1], ...
        'MarkerSize', 8);
end

if(strcmp(legendEntry, 'line'))
    % Makes the points not show up in the legend.
    set(get(get(h2, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
elseif(strcmp(legendEntry, 'point'))
    % Makes the lines not show up in the legend.
    set(get(get(h1, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
elseif(strcmp(legendEntry, 'none'))
    set(get(get(h1, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
    set(get(get(h2, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
end

if(~strcmp(location, 'None'))
    if(parSetColor(1) == 64/255 && parSetColor(2) == 64/255 && parSetColor(3) == 64/255)
        legend('Joint RSE', 'cursorGoal', 'Random Walk', 'Location', location);
    elseif(parSetColor(1) == 0 && parSetColor(2) == 189/255 && parSetColor(3) == 57/255)
        legend('Joint RSE', 'Lockstep RSE/RSE', 'Location', location);
    elseif(parSetColor(1) == 1 && parSetColor(2) == 133/255 && parSetColor(3) == 0)
        legend('Lockstep RSE/RSE', 'Lockstep RSE/RW', 'Location', location);
    end
end
