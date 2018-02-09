function [F0,JF,GJF,nit] = GOPT(J,GJ,F0,epsil,nitmax)
% Algorithme de descente du gradient à pas optimal avec intégration du
% calcul optimal par des formules analytiques pour les fonctions étudiées
%% Paramètres et déclarations utiles
nit = 0;
U = direct(F0);
Grad = adjoint(U);
error = norm(Grad);

%temporaire
x = linspace(0,1,1000);
solex = pi^2.*sin(pi.*x)';

%% On itère le processus jusqu'à se rapprocher de la solution
while (error > epsil) && (nit < nitmax)

    
    dk = -Grad;
    pas = pasopt(F0,dk);
    F0 = F0 + pas .* dk;
    
    U = direct(F0);
    Grad = adjoint(U);
    error = norm(Grad);
    nit = nit +1;
end
JF = J(F0);
GJF = GJ(F0);


end

