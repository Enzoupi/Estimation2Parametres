function [U] = direct(F)
%Fonction qui résout le problème (1) du TP2 avec la source F
n = length(F) + 1;
h = 1/n;
A = spmata(n-1);
U = h^2.*(A\F);
end

