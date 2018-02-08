function [A] = spmata(n)
%Fonction qui crée la matrice A comme définie dans le TP 2
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);
end

