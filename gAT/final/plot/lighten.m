function newcolor = lighten(color, amt)

% amt is
% 0 (default) for small (shading)
% 1 for big (gray out)

color = string_to_color(color);

if nargin < 2
   amt = 0;
end

if amt == 0
   p = 0.6;
elseif amt == 1
   p = 0.3;
end

newcolor = color*p + [1 1 1]*(1-p);

end
