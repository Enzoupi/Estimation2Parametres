function [F,JF,GJF,nit] = GCDYOPT(J,GJ,F0,epsil,nitmax)
%% Initialisation de Paramètres divers
nit = 0;            % Nombre d'itérations

%% Première étape : gradient à pas optimal normal
U = direct(F0);
Grad = adjoint(U);
error = norm(Grad);

if (error > epsil)
    dk = -Grad;
    % Calcul du pas optimal via la fonction dédiée
    pas = pasopt(F0,dk);
    % Mise à jour du premier point
    F = F0 + pas .* dk;
    U = direct(F);
    Grad_m1 = Grad;
    Grad = adjoint(U);
    error = norm(Grad);
    nit = nit +1;
    dkm1 = dk;
end

%% Gradient conjugué de Dai Yuan
while (error > epsil) && (nit < nitmax)
    % Calcul de Bkm1
    Bkm1 = sum( Grad .* Grad ) / ( sum(dkm1 .* (Grad-Grad_m1)) ) ;
    % Calcul du prochain pas de descente et MaJ du pas précédent
    dk = -Grad + Bkm1.*dkm1;
    % Calcul du pas optimal
    pas = pasopt(F,dk);
    % Calcul du nouveau F et mise à jour
    F = F + pas .* dk;
    U = direct(F);
    Grad_m1 = Grad;
    Grad = adjoint(U);
    % MàJ de la direction de la descente
    dkm1 = dk;
    error = max(abs(x-solex));
    nit = nit +1;
end
JF = J(F);
GJF = GJ(F);

