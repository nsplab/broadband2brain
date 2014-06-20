function thres = get_threshold(data_set, elec, id)
% Hard-coded eyeballed threshold values for analog thresholding

%if nargin < 3 || strcmp(id, 'ih')
    if data_set == 1
        if elec == 3
            thres = 100;
        elseif elec == 8
            thres = 75;
        elseif elec == 9
            thres = 200;
        elseif elec == 10
            thres = 110;
        elseif elec == 11
            thres = 60;
        elseif elec == 12
            thres = 75;
        elseif elec == 13
            thres = 35;
        elseif elec == 14
            thres = 100;
        elseif elec == 15
            thres = 80;
        elseif elec == 16
            thres = 125;
        end
    end
%end

end
