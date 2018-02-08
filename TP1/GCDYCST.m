function [x,Jx,GJx,nit] = GCDYCST(J,GJ,x0,pas,epsil,nitmax,findic)
% Initialisations
error = 10000;
nit = 0;

% Premiere etape Gradient à Pas Constant
x = x0 + pas .* (-GJ(x0,findic));
xm1 = x0;
dkm1 = -GJ(x0,findic);
nit = nit +1;
error = abs(norm(GJ(x,findic)));

%Gradient conjugue de Dai Yuan
while (error > epsil) && (nit < nitmax)
   Bkm1 = sum(GJ(x,findic).*GJ(x,findic))/(sum(dkm1.*(GJ(x,findic)-GJ(xm1,findic))));
   dk = -GJ(x,findic) + Bkm1.*dkm1;
   xm1 = x;
   x = x + pas .* dk;
   dkm1 = dk;
   error = abs(norm(GJ(x,findic)));
   nit = nit +1;
end
Jx = J(x,findic);
GJx = GJ(x,findic);