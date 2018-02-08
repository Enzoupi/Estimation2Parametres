clear variables
close all
%% Zone de commentaire
% Programme principal pour le TP1 d'estimation de paramètres

% remarque : question 1 regarder aussi la convergence et non pas que le
% nombre d'iterations

%Remarque : fonction de Rosenvrock : le gradient à pas constant va
%fonctionner pour des pas assez petits

%% Cas test avec un vecteur 
pas = 0.05;
epsil = 0.01;
nitmax = 1000;
findic = 1;
x0 = [2,2,2,2,5];
[x,Jx,GJx,nit] = GCST(@J,@GJ,x0,pas,epsil,nitmax,findic);

%% Question 2
figure
hold on
for ndim = 1
    % Solution exacte
    sol = 1:length(ndim);
    findic=1;
    x0 = zeros(1,ndim);
    epsil = 1e-7;
    pas = 0.1:0.1 :5;
    nitmax = 5000;
    nitconv=zeros(1,length(pas));
    res=zeros(1,length(pas));
    for i = 1:length(pas)
        [res(i),~,~,nitconv(i)]=GCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    end
    plot(pas,nitconv)
end
title('Convergence en fonction du pas')
xlabel('Taille du pas')
ylabel('Nombre diterations avant convergence')
legend('n=1')
hold off

figure
loglog(pas,abs(sol-res)/sol)
title('Ecart relatif en fonction du pas')
xlabel('Taille du pas')
ylabel('Ecart relatif')

%% Cas Test avec Dai Yuan
pas = 0.05;
epsil = 1e-7;
nitmax = 1000;
findic = 1;
x0 = [4,4,4];
[x,Jx,GJx,nit] = GCDYCST(@J,@GJ,x0,pas,epsil,nitmax,findic);

%% Question 3 : Convergence avec Dai Yuan
findic = 1;
pas = 0.1:0.1 :5;
err = zeros(2,length(pas));
res = zeros(1,length(pas));
nitmax = 5000;
figure
hold on
for ndim = 1:2
    solex = 1:ndim;
    findic=1;
    x0 = zeros(1,ndim);
    epsil = 1e-7;
    nitconv=zeros(1,length(pas));
    for i = 1:length(pas)
        [res,~,~,nitconv(i)]=GCDYCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    end
    err(ndim,i) = max(abs(solex-res)/solex);
    plot(pas,nitconv)
end
title('Convergence Pour GC de Dai Yuan')
xlabel('Taille du pas')
ylabel('Nombre diterations avant convergence')
legend('n=1','n=2')
hold off

%% Question 4 :
findic = 2;
x0 = [0,0];
solex = 1:length(x0);
pas = 0.1:0.1 :5;
errgcdy = zeros(1,length(pas));
errgcst = zeros(1,length(pas));
for i=1:length(pas)
    [resgcdy,~,~,nitgcdy(i)]=GCDYCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    [resgcst,~,~,nitg(i)]=GCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    errgcdy(i)=max(abs(solex-resgcdy)/solex);
    errgcst(i)=max(abs(solex-resgcst)/solex); 
end
figure
hold on
plot(pas,nitgcdy)
plot(pas,nitg)
title('Comparaison des convergences pour deux versions de la descente du gradient')
xlabel('Pas')
ylabel('Nombre d iterations avant convergence')
legend('Gradient conjugué (DY)','Gradient à pas constant')
figure
hold on
plot(pas,errgcdy)
plot(pas,errgcst)
ylim([0 1e-7])
title('Ecart relatif maximum')
xlabel('Pas')
ylabel('Ecart relatif')
legend('Gradient conjugué (DY)','Gradient à pas constant')

%% Question 5
findic = 3; % On choisit la fonction oscillante
x0 = 0:0.01:4;
pas = 0.01;
epsil = 1e-7;
nitmax = 5000;
% Détermination des minimas
res = zeros(1,length(x0));
for i=1:length(x0)
    [res(i),~,~,nitg(i)]=GCST(@J,@GJ,x0(i),pas,epsil,nitmax,findic);
end

figure
hold on
plot(x0,J(x0,findic))
plot(x0,res)
title('Minimum obtenu en fonction du point initial')
xlabel('x')
ylabel('y')
legend('Fonction oscillante','Minima obtenu')

%% Exercice 2 : Gradient à pas optimal
disp('---------------------------------------------------------------')
disp('Exercice 2')
disp('---------------------------------------------------------------')
% ==> Question 1 :

% Validation des descente de gradient à pas optimal : on sait que pour la
% fonction f on converge en une itération exactement et on retrouve ce
% comportement :
epsil = 1e-7;
nitmax = 1000;
findic = 1;
x0 = [5 3 64 15 -135];
disp('Avec le Gradient à pas optimal')
[x,Jx,GJx,nit,nbeval] = GOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Point de départ au hasard    : ' num2str(x0)])
disp(['Minimum obtenu               : ' num2str(x)])
disp(['En exactement ' num2str(nit) ' itération(s)'])
disp('Avec le Gradient Conjugué de Dai Yuan à pas optimal')
[x,Jx,GJx,nit,nbeval] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Point de départ au hasard    : ' num2str(x0)])
disp(['Minimum obtenu               : ' num2str(x)])
disp(['En exactement ' num2str(nit) ' itération(s)'])

%% Comparaison avec GCST pour n=2 et pas = 0.5
epsil = 1e-7;
nitmax = 1000;
findic = 1;
x0 = [15 -15];
pas = 0.5;
[resgopt,Jx,GJx,nit_gopt,nbeval] = GOPT(@J,@GJ,x0,epsil,nitmax,findic);
[resgcst,~,~,nit_gcst]=GCST(@J,@GJ,x0,pas,epsil,nitmax,findic);
[resgcdyopt,~,~,nit_gcdyopt,nbeval] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Comparaison sur la fonction g pour n = 2 avec une tolerance de ' num2str(epsil)])
disp(['Nombre d itérations pour GCST : ' num2str(nit_gcst)])
disp(['Nombre d itérations pour GOPT : ' num2str(nit_gopt)])
disp(['Nombre d itérations pour GCDYOPT : ' num2str(nit_gcdyopt)])


%% ==> Question 2 :
disp('Exercice 2 Question 2:')
% Pour f
%nn = [1,2,5,10,100,1000,10000];
nn = [1,2,3,4,5,6,7,8,9,10,12,15,17,20,25,50,75,100,250,500,750,1000,10000];
epsil = 1e-7;
nitmax = 10000;
findic =2;
err_gopt = zeros(size(nn));
err_gcdyopt = zeros(size(nn));
nit_gopt = zeros(size(nn));
nit_gcdyopt = zeros(size(nn));
nbeval_gopt = zeros(size(nn));
nbeval_gcdyopt = zeros(size(nn));
for i=1:length(nn)
    n=nn(i);
    x0 = zeros(1,n);
    solex = 1:n;
    [resgopt,~,~,nit_gopt(i),nbeval_gopt(i)]=GOPT(@J,@GJ,x0,epsil,nitmax,findic);
    [resgcdyopt,~,~,nit_gcdyopt(i),nbeval_gcdyopt(i)] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
    err_gopt(i) = max(abs(resgopt-solex));
    err_gcdyopt(i) = max(abs(resgcdyopt-solex));
end

err_gcdyopt(err_gcdyopt==0)=1e-17;
err_gopt(err_gopt==0)=1e-17;
figure
subplot(2,1,1)
hold on
plot(log10(nn),log10(err_gcdyopt),'-+')
plot(log10(nn),log10(err_gopt),'-+')
plot(log10(nn),ones(1,length(nn)).*-7)
xlabel('log10 de la dimension de l espace')
ylabel('log10 de l erreur')
legend('Gradient Conjugué DY Optimal','Gradient Optimal','Précision demandée','Location','northwest')
subplot(2,1,2)
hold on
plot(log10(nn),log10(nit_gcdyopt),'-+')
plot(log10(nn),log10(nit_gopt),'-+')
xlabel('log10 de la dimension de l espace')
ylabel('log10 du nombre d itérations')
legend('Gradient Conjugué DY Optimal','Gradient Optimal')



figure
hold on
plot(log10(nn),log10(nbeval_gcdyopt))
plot(log10(nn),log10(nbeval_gopt))
legend('Gradient Conjugué DY Optimal','Gradient Optimal','Location','Southeast')
xlabel('log10(n)')
ylabel('log10 nombre d évaluation de J et GJ')

%% représentation des trajets
epsil = 1e-7;
nitmax = 1000;
findic = 2;
x0 = [2.5 -2.5];
[resgopt,Jx,GJx,nit_gopt,nbeval] = GOPT_allx(@J,@GJ,x0,epsil,nitmax,findic);
[resgcdyopt,Jx,GJx,nit,nbeval] = GCDYOPT_allx(@J,@GJ,x0,epsil,nitmax,findic);


[X,Y]=meshgrid(-3:0.05:3);
Z=(X-1).^2+1/2.*(Y-2).^2
figure
hold on
plot3(resgopt(:,1),resgopt(:,2),(resgopt(:,1)-1).^2+1/2.*(resgopt(:,2)-2).^2,'-or')
plot3(resgcdyopt(:,1),resgcdyopt(:,2),(resgcdyopt(:,1)-1).^2+1/2.*(resgcdyopt(:,2)-2).^2,'-oy')
surf(X,Y,Z, 'edgecolor','none')
legend('Gradient optimal','Gradient Conjugué DY Optimal')
xlabel('x')
ylabel('y')

%% fonctions implémentées dans matlab
% fonctions transparentes pour définir au bon format f,g,s et h
func_f = @(x) J(x,1);
func_g = @(x) J(x,2);
func_s = @(x) J(x,3);
func_rosen = @(x) J(x,4);
options = optimset('TolX',1e-7,'Display','iter','MaxIter',50000,'MaxFunEvals',50000);
[xg,~,~,output_g] = fminsearch(func_g,[15,-15],options)
[xf,~,~,output_f] = fminsearch(func_f,[15,-15],options)
[xs,~,~,output_s] = fminsearch(func_s,[15],options)
[xrosen,~,~,output_rosen] = fminsearch(func_rosen,[15,35],options)
% Affichage des résultats
disp('--------------------------------------------------------------------')
disp('Exercice 3')
disp('--------------------------------------------------------------------')
disp('Résultats de minimisation via fminsearch sur f')
disp(['Minimum : ' num2str(xf)])
disp(['Erreur  : ' num2str(max(abs(xf-[1 2])))])
disp(['Nombre d itérations : ' num2str(output_f.iterations)])
disp('Résultats de minimisation via fminsearch sur g')
disp(['Minimum : ' num2str(xg)])
disp(['Erreur  : ' num2str(max(abs(xg-[1 2])))])
disp(['Nombre d itérations : ' num2str(output_g.iterations)])
disp('Résultats de minimisation via fminsearch sur s')
disp(['Minimum : ' num2str(xs)])
disp(['Erreur  : ' num2str(max(abs(xs-1)))])
disp(['Nombre d itérations : ' num2str(output_s.iterations)])
% 
% 
% % Test avec n=100
% x0 = zeros(1,10);
% [xg_big,~,~,output_g_big] = fminsearch(func_f,x0,options)
%% Minimisation de la fonction de Rosenbrock
disp('Minimisation de la fonction de Rosenbrock')
% Paramètres
x0=[0.8 1.5];
epsil = 1e-7;
nitmax = 100000;
findic = 4;
pas = 1e-4;
% calculs de minimisation
%[resgopt,~,~,nit_gopt,nbeval]   =   GOPT(     @J,@GJ,x0,      epsil,nitmax,findic);
%[resgcdyopt,~,~,nit,nbeval]     =   GCDYOPT(  @J,@GJ,x0,      epsil,nitmax,findic);
[resgcst,~,~,nit_gcst]             =   GCST(          @J,@GJ,x0,pas,  epsil,nitmax,findic);
[resgcdycst,~,~,nit_gcdycst]       =   GCDYCST(       @J,@GJ,x0,pas,  epsil,nitmax,findic);

findic = 4;
x0 = [0.8 1.5];
solex = [1,1];
pas = [1e-5 5e-5 1e-4 5e-4 1e-3 5e-3 1e-2 1e-1];
errgcdy = zeros(1,length(pas));
errgcst = zeros(1,length(pas));
nitgcdy = zeros(1,length(pas));
nitg = zeros(1,length(pas));
for i=1:length(pas)
    [resgcdy,~,~,nitgcdy(i)]=GCDYCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    [resgcst,~,~,nitg(i)]=GCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    errgcdy(i)=max(abs(solex-resgcdy)/solex);
    errgcst(i)=max(abs(solex-resgcst)/solex); 
end
figure
hold on
plot(log10(pas),log10(nitgcdy),'-+')
plot(log10(pas),log10(nitg),'-+')
title('Comparaison des convergences pour deux versions de la descente du gradient')
xlabel('log10 du pas')
ylabel('log10 du nombre d iterations avant convergence')
legend('Gradient conjugué (DY)','Gradient à pas constant')
figure
hold on
plot(log10(pas),log10(errgcdy),'-+')
plot(log10(pas),log10(errgcst),'-+')
title('Ecart relatif maximum')
xlabel('log10 du pas')
ylabel('log10 de l écart relatif')
legend('Gradient conjugué (DY)','Gradient à pas constant')