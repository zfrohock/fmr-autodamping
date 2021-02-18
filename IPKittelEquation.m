%% The in plane Kittel equation
function [y] = IPKittelEquation(x,gamma,Ms)
% The in plane Kittel equation. Y is FMR field, x is FMR frequency
y = (-Ms + sqrt((Ms).^2 + 4.*x.^2./gamma^2))./2;
end