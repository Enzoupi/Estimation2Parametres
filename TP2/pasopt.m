function [rho] = pasopt(F,D)
% PASOPT fonction qui calcule le pas optimal en F dans la direction D

% 1) Calcul de W et U
W = differentielle(D);
U = direct(F);

% 2) Appel a la variable globale Uobs (Danger potentiel !)
global Uobs

% 3) Calcul du pas optimal
% Norme L2 : Il apparaît un h en haut et en bas, donc ils se simplifient
Norme_L2_W = sum(W.*W) ;
rho = - sum( (U-Uobs) .* W) / Norme_L2_W;

%|----------------------------------------------------------------
%| Vérifier que le produit scalaire dans L2 prend lui aussi un h (a priori
%| oui parce qu'on intègre !)
%|----------------------------------------------------------------

end

