function [F,JF,GJF,nit] = GCST(J,GJ,F0,pas,epsil,nitmax)
%% Calcul de lerreur au point initial
nit = 0;
Grad = GJ(F0);
error = norm(Grad);
F = F0;
while (error > epsil) && (nit < nitmax)
    dk = -Grad;
    F = F + pas .* dk;
    Grad = GJ(F);
    error = norm(Grad);
    nit = nit +1;
end
JF = J(F);
GJF = GJ(F);

