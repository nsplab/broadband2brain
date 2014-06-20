function newColor = lightenColor(color)

    if(color(1) == 1 && color(2) == 35/255 && color(3) == 0)
        newColor = [1 90/255 64/255];
    elseif(color(1) == 7/255 && color(2) == 118/255 && color(3) == 160/255)
        newColor = [59/255 165/255 211/255];
    elseif(color(1) == 64/255 && color(2) == 64/255 && color(3) == 64/255)
        newColor = [128/255 128/255 128/255];
    elseif(color(1) == 0 && color(2) == 189/255 && color(3) == 57/255)
        newColor = [56/255 223/255 113/255];
    elseif(color(1) == 1 && color(2) == 133/255 && color(3) == 0)
        newColor = [1 163/255 64/255];
    else
        newColor = color;
    end
    
end
