%% Equation used for Gaussian fitting
function [y] = GaussianDerivativeFunction(b,c,height,Hpp,Hfmr,x)
x = x';
y = zeros(length(x));
y = b + c.*x + -1*((Hpp/2)*exp(.5))*(height*(-x+Hfmr).*exp(-((-x+Hfmr).^2)/(2*(Hpp/2)^2)))/((Hpp/2)^2);
y = y';
end