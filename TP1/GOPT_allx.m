function [x,Jx,GJx,nit,counter] = GOPT_allx(J,GJ,x0,epsil,nitmax,findic)
% Algorithme de descente du gradient à pas optimal avec intégration du
% calcul optimal par des formules analytiques pour les fonctions étudiées
%% Paramètres et déclarations utiles
error = 10000;
nit = 0;
i = 1:length(x0);
counter = 0; % compteur d'évaluation de la fonction et de son gradient
% solution exacte pour calcul de l'erreur
if findic == 1 || findic == 2
    solex = i;
end
% Vecteur pour conserver tous les x en sortie de l'algo
x=zeros(nitmax,length(x0));
x(1,:) = x0;
%% On itère le processus jusqu'à se rapprocher de la solution
while (error > epsil) && (nit < nitmax)
    dk = -GJ(x0,findic); counter = counter +1;
    if findic == 1
        b = 2*sum( (x0-i) .* dk);
        c = sum(dk .* dk);
        pas = -b/(2*c);
    elseif findic == 2
        b = 2*sum( 1./i .* (x0-i) .* dk);
        c = sum(1./i .* dk .* dk);
        pas = -b/(2*c);
    else
        disp('---------------------------------------------------------------')
        disp('La fonction demandée nest pas encore implémentée !!')
        disp('---------------------------------------------------------------')
        return
    end
    x0 = x0 + pas .* dk;
    error = max(abs(x0-solex));
    nit = nit +1;
    x(nit+1,:)=x0;
end
Jx = J(x0,findic);
GJx = GJ(x0,findic);
counter = counter +2;

x = x(1:nit+1,:);
end



