function [x,Jx,GJx,nit,counter] = GCDYOPT(J,GJ,x0,epsil,nitmax,findic)
%% Paramètres divers
error = 10000;      % Initialisation d'une valeur élévée pour l'erreur
nit = 0;            % Nombre d'itérations
i = 1:length(x0);   % Comodité pour réaliser les calculs
counter = 0;        % Compteur dévaluation de la fonction et de son gradient


%% Première étape : gradient à pas optimal normal
% Initialisation des premiers calculs de pas optimaux
dk = -GJ(x0,findic);
counter = counter +1;
if findic == 1
    b = 2*sum( (x0-i) .* dk);
    c = sum(dk .* dk);
    pas = -b/(2*c);
    solex = i;
elseif findic == 2
    b = 2*sum( 1./i .* (x0-i) .* dk);
    c = sum(1./i .* dk .* dk);
    pas = -b/(2*c);
    solex = i;
elseif findic == 4
    solex = [1 1];
    % FAUX !!!
    disp('Attention : le calcul du pas optimal n est pas bon !')
    pas = (1-x0(1))/(dk(1));
else
    disp('---------------------------------------------------------------')
    disp('La fonction demandée nest pas encore implémentée !!')
    disp('---------------------------------------------------------------')
    return
end
x = x0 + pas .* dk;
error = max(abs(x-solex));
nit = nit +1;

dkm1 = dk;
xm1 = x0;
GJx = -dk;

%% Gradient conjugué de Dai Yuan
while (error > epsil) && (nit < nitmax) %&& (norm(GJx) > epsil)
    GJx = GJ(x,findic);
    counter = counter +1;
    Bkm1 = sum( GJx .* GJx ) / ( sum(dkm1 .* (GJx-GJ(xm1,findic)) ) );
    dk = -GJx + Bkm1.*dkm1;
    % Calcul du pas optimal
    if findic == 1
        b = 2 * sum( (x-i) .* dk );
        c = sum( dk .* dk );
        pas = -b/(2*c);
    elseif findic == 2
        b = sum( 1./i .* 2 .* (x-i) .* dk);
        c = sum(1./i .* dk .* dk);
        pas = -b/(2*c);
    elseif findic == 4
        % FAUX !!!
        disp('Attention : le calcul du pas optimal n est pas bon !')
        pas = (1-x0(1))/(dk(1));
    end
    xm1 = x;
    x = x + pas .* dk;
    
    dkm1 = dk;
    error = max(abs(x-solex));
    nit = nit +1;
end
Jx = J(x,findic);
GJx = GJ(x,findic);
counter = counter +2;

