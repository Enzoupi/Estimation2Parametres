function [allx,Jx,GJx,nit] = GCDYCST_allx(J,GJ,x0,pas,epsil,nitmax)
% Fonction qui procède à la minimisation d'une fonction par l'algorithme du
% gradient. C'est une variante qui possède en sortie toutes les itérations
% successives de l'algortihme. Pour un gain en performance, utiliser plutôt
% GCDYCST.

% Initialisations diverses
error = 10000;
nit = 0;

% Matrice pour retenir toutes les itérations
allx(:,1) = x0;

% Premiere etape Gradient à Pas Constant
x = x0 + pas .* (-GJ(x0));
xm1 = x0;
dkm1 = -GJ(x0);
nit = nit +1;
allx(:,nit+1) = x;
error = abs(norm(GJ(x)));

%Gradient conjugue de Dai Yuan
while (error > epsil) && (nit < nitmax)
   Bkm1 = sum(GJ(x).*GJ(x))/(sum(dkm1.*(GJ(x)-GJ(xm1))));
   dk = -GJ(x) + Bkm1.*dkm1;
   xm1 = x;
   x = x + pas .* dk;
   dkm1 = dk;
   error = abs(norm(GJ(x)));
   nit = nit +1;
   allx(:,nit+1) = x;
end
Jx = J(x);
GJx = GJ(x);