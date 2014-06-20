function h = scatterBeautifulData(xData, yData, parSetColor, location)

if(nargin < 4)
    location = 'Best';
end

hold on;
h = scatter(xData, yData, 'MarkerEdgeColor', parSetColor, 'LineWidth', 2, ...
    'Marker', 'o', 'MarkerFaceColor', lightenColor(parSetColor), ...
    'SizeData', 64);

if(~strcmp(location, 'None'))
    if(parSetColor(1) == 7/255 && parSetColor(2) == 118/255 && parSetColor(3) == 160/255)
        legend('Joint RSE', 'cursorGoal', 'Location', location);
    elseif(parSetColor(1) == 64/255 && parSetColor(2) == 64/255 && parSetColor(3) == 64/255)
        legend('Joint RSE', 'cursorGoal', 'Random Walk', 'Location', location);
    elseif(parSetColor(1) == 0 && parSetColor(2) == 189/255 && parSetColor(3) == 57/255)
        legend('Joint RSE', 'Lockstep RSE/RSE', 'Location', location);
    elseif(parSetColor(1) == 1 && parSetColor(2) == 133/255 && parSetColor(3) == 0)
        legend('Lockstep RSE/RSE', 'Lockstep RSE/RW', 'Location', location);
    end
end
