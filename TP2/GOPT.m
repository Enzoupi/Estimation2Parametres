function [F0,JF,GJF,nit,allpasopt] = GOPT(J,GJ,F0,epsil,nitmax)
% Algorithme de descente du gradient à pas optimal avec intégration du
% calcul optimal par des formules analytiques pour les fonctions étudiées
%% Paramètres et déclarations utiles
nit = 0;
Grad = GJ(F0);
error = norm(Grad);

%% On itère le processus jusqu'à se rapprocher de la solution
while (error > epsil) && (nit < nitmax)
    dk = -Grad;
    pas = pasopt(F0,dk);
    F0 = F0 + pas .* dk;
    Grad = GJ(F0);
    error = norm(Grad);
    nit = nit +1;
    allpasopt(nit)=pas;
end
JF = J(F0);
GJF = GJ(F0);


end

