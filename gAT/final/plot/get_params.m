function [style, width, markersize] = get_params(dots, thick)

if dots == 0;
   style = '-';
elseif dots == 1
   style = '--';
end

if thick == 0
   width = 5;
   markersize = 12;
elseif thick == 1
   width = 14;
   markersize = 30;
elseif thick == 2
    width = 2;
    markersize = 0;
end


end