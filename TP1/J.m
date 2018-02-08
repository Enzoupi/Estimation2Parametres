function [resJ] = J(x,findic)
%f(x)
if (findic == 1)
    i=1:length(x);
    resJ = sum((x-i).*(x-i));
% g(x)
elseif (findic == 2)
    i=1:length(x);
    resJ = sum(1./i.*(x-i).*(x-i));
% Fonction oscillante s
elseif (findic == 3)
    resJ=(x-1).^2+sin(4*pi*(x-9/8));
% Fonction de Rosenbrock
elseif (findic == 4)
    if length(x) ~= 2
        disp('X does not have the proper size for Rosenbrock')
    end
    resJ=(1-x(1)).^2 + 100.*(x(2)-x(1).^2).^2;
end