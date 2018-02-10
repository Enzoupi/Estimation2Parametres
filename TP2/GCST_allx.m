function [allx,Jx,GJx,nit] = GCST_allx(J,GJ,x0,pas,epsil,nitmax)
error = max(abs(GJ(x0)));
x = x0;
nit = 0;
allx(:,1)=x;
while (error > epsil) && (nit < nitmax)
   x = x - pas .* GJ(x);
   nit = nit +1
   allx(:,nit+1)=x;
   Jx = J(x);
   GJx = GJ(x);
   error = max(abs(GJx));
end



