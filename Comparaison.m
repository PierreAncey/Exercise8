%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=16;lw=2;
fichier = 'output';
data = load([fichier,'_observables.out']);
t = data(:,1);
xmoy = data(:,6);
pmoy = data(:,8);

x = 64.571596*sin(0.004*t);
p = 0.258286*cos(0.004*t);

%% Figures %%
%%%%%%%%%%%%%

figure
plot(t,xmoy,'-b','linewidth',lw)
hold on
plot(t,x,'--r','linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('<x>')
legend('<x_{numerical}>','x_{classical}', 'Location', 'se')

figure
plot(t,pmoy,'-b','linewidth',lw)
hold on
plot(t,p,'--r','linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('<p>')
legend('<p_{numerical}>','p_{classical}', 'Location', 'se')

xbis = x;
xmoybis = xmoy;
pmoybis = pmoy;
pbis = p;
tbis = t;
for i=size(x):1
    if (abs(x(i)) < 1) || (abs(p(i)) < 0.004)
        xbis(i) =[];
        xmoybis(i) =[];
        pmoybis(i) =[];
        pbis(i) =[];
        tbis(i) =[];
    end
end
figure
plot(t,abs((xmoy-x)./xmoy),'-b','linewidth',lw)
hold on
plot(t,abs((pmoy-p)./pmoy),'--r','linewidth',lw)
set(gca,'fontsize',fs)
grid on
xlabel('t')
ylabel('Relative error')
leg = legend('$\left|(\frac{<x>-x_{classical}}{<x>})\right|$','$\left|(\frac{<p>-p_{classical}}{<p>})\right|$','Interpreter','latex','Location', 'nw');

