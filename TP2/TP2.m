clear all
close all
%% Zone de Commentaire
% Script pour le TP2 d'estimations de paramètres

%% Tests de convergence du problème direct

%Test dans un cas simple
n = 150;
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

%% Exercice 3 Gradient pas constant
% ==> Question 2) Faire des tests pour plusieurs valeurs du pas
n = 150;
h = 1/n;
x=linspace(h,1-h,n-1);
Uobs = sin(pi.*x)';
F = zeros(length(x),1);
solex = pi^2.*sin(pi.*x)';
epsil = 1e-7;
nitmax = 10000;
pas = [0.5 2 4 8 16 32 64 128 256 512];
err = zeros(1,length(pas));
niter = zeros(1,length(pas));
for i=1:length(pas)
    [xmin,Jxmin,GJxmin,nit] = GCST(@J,@GJ,F,pas(i),epsil,nitmax);
    err(i) = max(abs(xmin-solex));
    niter(i)= nit;
    if pas(i) == 64
        xmin_64=xmin;
    end
    
end
figure
hold on
plot(x,xmin_64,'-+')
plot(x,solex)
legend('Solution approchée','Solution exacte')

figure
subplot(2,1,1)
plot(log10(pas),log10(niter),'-+')
xlabel('log10 du pas')
ylabel('log10 du nombre d itérations')
subplot(2,1,2)
plot(log10(pas),log10(err),'-+')
xlabel('log10 du pas')
ylabel('log10 de l erreur')

%% ====> Question 2 : Pour pas = 64 vérifier les formules analytiques
% dérivées en TD
pas = 64;
%pas = pi^2;
nn = [1e2 1e3 1e4];
epsil = 1e-7;
nitmax = 10000;
nitmax_effect = 0;
err = zeros(length(nn),nitmax);
for i=1:length(nn)
    n = nn(i);
    h = 1/n;
    x=linspace(h,1-h,n-1);
    Uobs = sin(pi.*x)';
    F = zeros(length(x),1);
    solex = pi^2.*sin(pi.*x)';
    
    [xmin,Jxmin,GJxmin,nit] = GCST_allx(@J,@GJ,F,pas,epsil,nitmax);
    nitmax_effect = max(nitmax_effect,nit);
    err(i,1:(nit+1))=max(abs(xmin-solex));
end
err=err(:,1:nitmax_effect+1);

%Calcul de l'ordre de convergence
coefs = polyfit(1:(nitmax_effect-3),log10(err(end,2:end-3)),1);

figure
hold on
plot(0:nitmax_effect,log10(err),'-+')
plot(0:nitmax_effect,coefs(1).*(0:nitmax_effect)+coefs(2),'--')
legend('n=100','n=1000','n=10000',[num2str(coefs(1)) '* iter + ' num2str(coefs(2))])
xlabel('Itération')
ylabel('log10 de l''erreur')
title('Erreur en fonction de l''itération')

figure
hold on
plot(x,xmin(:,[1,2,3,4,10,end]))
plot(x,solex,'--k')
legend('iteration 0','iteration 1','iteration 2','iteration 3','iteration 9',['iteration ' num2str(nit+1)],'solution exacte')
xlabel('x')
ylabel('y')

%% Exercice 3 : Gradient conjugué à pas constant
n = 100;
h = 1/n;
x=linspace(h,1-h,n-1);
Uobs = sin(pi.*x)';
F = zeros(length(x),1);
solex = pi^2.*sin(pi.*x)';
epsil = 1e-7;
nitmax = 10000;
pas = [0.5 2 4 8 16 32 64 128 256 512];
err = zeros(1,length(pas));
niter = zeros(1,length(pas));
for i=1:length(pas)
    [xmin,Jxmin,GJxmin,nit] = GCDYCST(@J,@GJ,F,pas(i),epsil,nitmax);
    err(i) = max(abs(xmin-solex));
    niter(i)= nit;
    if pas(i) == 64
        xmin_64=xmin;
    end
    
end
figure
hold on
plot(x,xmin_64,'-+')
plot(x,solex)
legend('Solution approchée','Solution exacte')
title('GCDYCST')

figure
subplot(2,1,1)
plot(log10(pas),log10(niter),'-+')
xlabel('log10 du pas')
ylabel('log10 du nombre d itérations')
subplot(2,1,2)
plot(log10(pas),log10(err),'-+')
xlabel('log10 du pas')
ylabel('log10 de l erreur')


%% ====> Question 2 : Pour pas = 64 vérifier les formules analytiques
% dérivées en TD
pas = 64;
%pas = pi^4;
nn = [1e2 1e3 1e4];
epsil = 1e-7;
nitmax = 100;
nitmax_effect = 0;
err = zeros(length(nn),nitmax);
for i=1:length(nn)
    n = nn(i);
    h = 1/n;
    x=linspace(h,1-h,n-1);
    Uobs = sin(pi.*x)';
    F = zeros(length(x),1);
    solex = pi^2.*sin(pi.*x)';
    
    [xmin,Jxmin,GJxmin,nit] = GCDYCST_allx(@J,@GJ,F,pas,epsil,nitmax);
    nitmax_effect = max(nitmax_effect,nit);
    err(i,1:(nit+1))=max(abs(xmin-solex));
end
err=err(:,1:nitmax_effect+1);

figure
hold on
plot(0:nitmax_effect,log10(err)','-+')
legend('n=100','n=1000','n=10000')
xlabel('Itération')
ylabel('log10 de l''erreur')
title('Erreur en fonction de l''itération')

figure
hold on
plot(x,xmin(:,[1,2,end]))
plot(x,solex,'--k')
legend('Initialisation F(0)','iteration 1',['iteration ' num2str(nit)],'solution exacte')
xlabel('x')
ylabel('y')

%% Exercice 4 : Calcul de la différentielle
% Sur quel exemple peut-on valider ? C'est le même probleme que le problème
% direct de la première question. On peut faire le test avec la même
% fonction.
x = linspace(0,1,1e4);
D = pi^2*sin(pi.*x)';
solex = sin(pi.*x)';
W = differentielle(D);
figure
hold on
plot(x,solex)
plot(x,W,'--')
legend('solex','W')
title('Solution exacte et approchée superposées')

figure
semilogy(x,abs(W-solex))
title('log de l''erreur en fonction de la discrétisation en espace')
xlabel('Discrétisation en espace')
ylabel('log de l''erreur')

%% Exercice 5 Question 2)
disp('Exercice 5 question 2')
epsil=1e-7;
nitmax=1000;
nn=[100,250,500,750,1000,5000,10000];
err=zeros(size(nn));
for i=1:length(nn)
    n=nn(i);
    h = 1/n;
    x=linspace(h,1-h,n-1);
    solex = pi^2.*sin(pi.*x)';
    Uobs = sin(pi.*x)';
    F0=zeros(size(x))';

    [F,JF,GJF,nit]=GOPT(@J,@GJ,F0,epsil,nitmax);
    err(i) = max(abs(solex-F));
end
figure
plot(x,F)

%Calcul de la convergence
coefs = polyfit(log10(nn),log10(err),1)

figure
loglog(nn,err,'-+')
xlabel('Nombre de points de discrétisation de l''espace')
ylabel('maximum de l''erreur')
title(['Convergence de GOPT : ' num2str(-coefs(1))])

%% Exercice 5 Question 3)
disp('Exercice 5 question 3')
epsil=1e-7;
nitmax=1000;
nn=[5,10,20,50,75,100,250,500,750,1000,5000,10000];
err=zeros(size(nn));
figure
hold on
for i=1:length(nn)
    n=nn(i);
    h = 1/n;
    x=linspace(h,1-h,n-1);
    solex = pi^2.*sin(pi.*x)';
    Uobs = sin(pi.*x)';
    F0=zeros(size(x))';
    epsil=1e-7;
    nitmax=100;
    [F,JF,GJF,nit]=GCDYOPT(@J,@GJ,F0,epsil,nitmax);
    err(i) = max(abs(solex-F));
    plot(x,F)
end
xlabel('abscisse')
ylabel('ordonnées')
title('Résultats obtenus pour différents n')

%Calcul de la convergence
coefs = polyfit(log10(nn),log10(err),1)

figure
loglog(nn,err,'-+')
xlabel('Nombre de points de discrétisation de l''espace')
ylabel('maximum de l''erreur')
title(['Convergence de GCDYOPT : ' num2str(-coefs(1))])

