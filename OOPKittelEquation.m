%% The out of plane Kittel equation
function [y] = OOPKittelEquation(x,gamma,Ms)
% The out of plane Kittel equation. X is operational frequency, Y is FMR
% field
y = x./gamma + 4*pi*Ms;
end