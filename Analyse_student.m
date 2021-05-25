%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output';
data = load([fichier,'_observables.out']);
t = data(:,1);
Pgauche = data(:,2);
Pdroite = data(:,3);
Ptot = data(:,4);
E = data(:,5);
xmoy = data(:,6);
x2moy = data(:,7);
pmoy = data(:,8);
p2moy = data(:,9);
DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
Deltax = data(:,11);
Deltap = data(:,12);
data = load([fichier,'_potential.out']);
x = data(:,1);
V = data(:,2);
tpsi2 = load([fichier,'_psi2.out']);
[nt,nx1]=size(tpsi2);
psi2 = tpsi2(2:nt,2:nx1);

%% Figures %%
%%%%%%%%%%%%%
fs=16;lw=2; % fs est utilisé pour la taille de la police, lw pour l'épaisseur des lignes
figure
plot(t,Pgauche, t,Pdroite,t,Ptot,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('P')
legend('P_{x<xa}','P_{x>xb}','P_{tot}','Location','Best')
%figure('Name',['Analyse de ' fichier], 'unit', 'norm', 'pos', [1,.1,.8,.3])

%subplot(2,2,1)
figure
plot(x,V,'linewidth',lw)
set(gca,'fontsize',fs)
hold on
plot(x([1,end]),E(1)*ones(1,2),'--','linewidth',lw)
grid on
xlabel('x')
ylabel('V')
legend('V(x)','E','Location','Best')

figure
[X,T] = meshgrid(x,t);
pcolor(X,T,psi2); % voir aussi contour, contourf
set(gca,'fontsize',fs)
shading interp
colormap jet % il existe d'autres colormaps (help colormap)
c = colorbar;
xlabel('x [m]')
ylabel('t [s]')
ylabel(c,'|\psi|^2')

figure
plot(t,E,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('E')

figure
plot(t,Deltax.*Deltap,'linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('<\Delta x>(t) \cdot <\Delta p>(t)')
