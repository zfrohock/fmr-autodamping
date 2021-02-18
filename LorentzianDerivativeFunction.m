%% Equation defining the Lorentzian fit
function [y] = LorentzianDerivativeFunction(A,b,c,Hpp,Hfmr,theta,x)
% Function to which we fit the FMR power derivative. 
x = x';
y = zeros(length(x));
y = b + c.*x + (-16).*A.*cos(theta.*pi./180).*Hpp.*sqrt(3).*(x-Hfmr)./pi./(4.*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2).^2 ...
    +   4.*A.*sin(theta.*pi./180).*((-4).*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2)/pi/(4.*(x-Hfmr).^2+(Hpp.*sqrt(3)).^2).^2;
y = y';
end