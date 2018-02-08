function [xres,Jx,GJx,nit,counter] = GCDYOPT_allx(J,GJ,x0,epsil,nitmax,findic)
%% Paramètres divers
error = 10000;      % Initialisation d'une valeur élévée pour l'erreur
nit = 0;            % Nombre d'itérations
i = 1:length(x0);   % Comodité pour réaliser les calculs
counter = 0;        % Compteur dévaluation de la fonction et de son gradient

% Vecteur pour conserver tous les x en sortie de l'algo
xres=zeros(nitmax,length(x0));
xres(1,:) = x0;


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
else
    disp('---------------------------------------------------------------')
    disp('La fonction demandée nest pas encore implémentée !!')
    disp('---------------------------------------------------------------')
    return
end
x = x0 + pas .* dk;
error = max(abs(x-solex));
nit = nit +1;
xres(nit +1,:)=x;

dkm1 = dk;
xm1 = x0;

%% Gradient conjugué de Dai Yuan
while (error > epsil) && (nit < nitmax)
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
        b = 2 * sum(1./i .* (x0-i) .* dk);
        c = sum(1./i.* dk .* dk);
        pas = -b/(2*c);
    end
    xm1 = x;
    x = x + pas .* dk;
    dkm1 = dk;
    error = max(abs(x-solex));
    nit = nit +1;
    xres(nit +1,:)=x;
end
Jx = J(x,findic);
GJx = GJ(x,findic);
counter = counter +2;

xres = xres(1:nit+1,:);

end