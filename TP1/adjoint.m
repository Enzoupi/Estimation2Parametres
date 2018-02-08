function [V] = adjoint(U)
global Uobs
% Fonction qui résout le problème adjoint comme mentionné dans l'exerice 2
% du TP 2
n = length(U) + 1;
h = 1/n;
A = spmata(n-1);
V=(h^2)*(A\(U-Uobs));
end

