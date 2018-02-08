function [x,Jx,GJx,nit] = GCST(J,GJ,x0,pas,epsil,nitmax,findic)
error = 10000;
nit = 0;
while (error > epsil) && (nit < nitmax)
   error = abs(norm(GJ(x0,findic)));
   x = x0 - pas .* GJ(x0,findic);
   nit = nit +1;
   x0=x;
   Jx = J(x,findic);
   GJx = GJ(x,findic);
end

