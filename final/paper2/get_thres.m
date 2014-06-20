% Thresholds for AT (analog thresholding) and aFRI
function t = get_thres(data_set, channel)

if data_set == 1
    
    % Optimal threshold for T = 0.005 on aFRI (gAT)
    % NOTE: Computed for data_set 1 => tests on data_set 1 are "overfit"

    if channel == 3
        t = 108.9;
    elseif channel == 8
        t = 81.46;
    elseif channel == 9
        t = 116.35;
    elseif channel == 10
        t = 111.47;
    elseif channel == 11
        t = 54.41;
    elseif channel == 12
        t = 104.27;
    elseif channel == 13
        t = 58.39;
    elseif channel == 14
        t = 114.3;
    elseif channel == 15
        t = 77.93;
    elseif channel == 16
        t = 113.09;
    end

elseif data_set == 2

    if channel == 1
        t = 47.298;
    elseif channel == 2
        t = 59.7312;
    elseif channel == 4
        t = 46.3835;
    elseif channel == 5
        t = 55.4064;
    elseif channel == 12
        t = 50.2258;
    elseif channel == 13
        t = 63.5136;
    elseif channel == 14
        t = 49.498;
    end

end
