clear all
close all
%% Zone de Commentaire
% Script pour le TP2 d'estimations de paramètres

%% Tests de convergence du problème direct

%Test dans un cas simple
n = 150
h = 1/n;
x=linspace(h,1-h,n-1);
F = pi^2*sin(pi.*x)';
solex = sin(pi.*x)';
U=direct(F);
err=sum(abs(U-solex));
figure
hold on
plot(x,solex,'--')
plot(x,F)
plot(x,U,'+')
legend('solution exacte','F','solution calculée')

nn=[10,20,30,40,50,100,250,500,1e3,1e4];
h=zeros(size(n));
err=zeros(size(nn));
for i=1:length(nn)
    n = nn(i); 
    h(i) = 1/n;
    x=linspace(h(i),1-h(i),n-1);
    F = pi^2*sin(pi.*x)';
    solex = sin(pi.*x)';
    U=direct(F);
    err(i)=max(abs(U-solex));
end
figure
loglog(h,err)
coeff = polyfit(log10(h),log10(err),1);
xlabel('log10 de h')
ylabel('log10 de l erreur')
title(['Convergence du probleme direct : ' num2str(coeff(1))])
U_direct = U;

%% Exercice 2 : Convergence du problème adjoint
nn=[10,20,30,40,50,100,250,500,1e3,1e4];
h=zeros(size(n));
err=zeros(size(nn));
global Uobs
for i=1:length(nn)
    n = nn(i); 
    h(i) = 1/n;
    x=linspace(h(i),1-h(i),n-1);
    %---------------------------------------------------------------------
%     % Test avec le U du sujet : Fonctionnel !
%     U = (pi^2+1)*sin(pi.*x);
    %---------------------------------------------------------------------
    % Test avec le U de U=direct(F) On doit changer le F pour avoir le bon
    % resultat : Fonctionnel aussi
    F = pi^2*(pi^2+1)*sin(pi.*x)'; 
    U = direct(F);
    %---------------------------------------------------------------------
    Uobs = sin(pi.*x)';
    V=adjoint(U);
    err(i)=max(abs(V-Uobs));
end
figure
loglog(h,err)
coeff = polyfit(log10(h),log10(err),1);
xlabel('log10 de h')
ylabel('log10 de l erreur')
title(['Convergence du probleme adjoint : ' num2str(coeff(1))])

%% Exercice 3


