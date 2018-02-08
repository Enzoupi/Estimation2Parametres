function [W] = differentielle(D)
%differentielle : calcule F'(f).d avec d un vecteur colonne de taille (n-1)
% ou n est la dimension de l'espace
n = length(D) + 1;
h = 1/n;
A = spmata(n-1);
W=(h^2)*(A\D);
end

