% Convexify an ROC curve with data points X, Y
function [x y] = convexify(X, Y)

x = zeros(size(X));
y = zeros(size(Y));

% Reverse if necessary
if median(diff(Y)) > 0
    X = X(end:-1:1);
   	Y = Y(end:-1:1);
end

% Make x increasing
for i = 1 : length(X)
    x(i) = min(X(i:end));
end

% Make y decreasing
for i = 1 : length(X)
    y(i) = min(Y(1:i));
end

end