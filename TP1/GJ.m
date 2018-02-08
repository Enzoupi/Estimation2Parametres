function [resGJ] = GJ(x,findic)
if (findic == 2)
   i=1:length(x);
   resGJ = 2./i.*(x-i);
elseif (findic == 1)
   i=1:length(x);
   resGJ = 2.*(x-i);
% Gradient Fonction oscillante
elseif (findic == 3)
    resGJ=2*(x-1)+4*pi*cos(4*pi*(x-9/8));
% Gradient de la fonction de Rosenbrock
elseif (findic == 4)
    resGJ(1) = 400*x(1)^3-2*x(1)*(200*x(2)-1)-2;
    resGJ(2) = 100*2*(x(2)-x(1)^2);
end
