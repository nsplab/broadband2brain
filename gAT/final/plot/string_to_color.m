function newcolor = string_to_color(color)

newcolor = color;

if isstr(color)

   if strcmpi(color, 'real')
      newcolor = rgbconv('666666');
   elseif strcmpi(color, 'gibbs')
      newcolor = rgbconv('06799F');
   elseif strcmpi(color, 'iterml')
      newcolor = rgbconv('FF8300');
   elseif strcmpi(color, 'err')
      newcolor = rgbconv('F10026');
   elseif strcmpi(color, 'ga')
       newcolor = rgbconv('3BDA00');
   end

end

end