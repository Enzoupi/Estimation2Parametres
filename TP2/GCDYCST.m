function [x,Jx,GJx,nit] = GCDYCST(J,GJ,x0,pas,epsil,nitmax)
% Initialisations
error = 10000;
nit = 0;

% Premiere etape Gradient à Pas Constant
x = x0 + pas .* (-GJ(x0));
xm1 = x0;
dkm1 = -GJ(x0);
nit = nit +1;
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
end
Jx = J(x);
GJx = GJ(x);