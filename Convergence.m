% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %% (A MODIFIER SELON VOS BESOINS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = 2; fs = 16; ms = 10;

repertoire = ''; % Chemin d'acces au code compile
executable = 'Exercice8'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

nsimul = 5; % Nombre de simulations a faire

dt = [0.0625 0.125 0.25 0.5 1];


paramstr = 'dt'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul);

for i = 1:nsimul
    filename = [paramstr, '=', num2str(param(i))];
    output{i} = [filename, '.out'];
    c = sprintf('%s %s %s %s%s %s%s %s=%.15g output_observables=%s', 'set', 'path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=%', '&', repertoire, executable, repertoire, input, paramstr, param(i), output{i})
    system(c);
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

x = zeros(1,nsimul);
x2 = zeros(1,nsimul);
p = zeros(1,nsimul);
p2 = zeros(1,nsimul);
error_x = zeros(1,nsimul);
error_p = zeros(1,nsimul);
error_xbis = zeros(1,nsimul);
error_pbis = zeros(1,nsimul);

for i = 1:nsimul
    data = load(output{i});
    x_temp = data(:,6);
    x2_temp = data(:,7);
    p_temp = data(:,8);
    p2_temp = data(:,9);
    error_xbis_temp = data(:,11);
    error_pbis_temp = data(:,12);
    
    x(i) = x_temp(end);
    x2(i) = x2_temp(end);
    p(i) = p_temp(end);
    p2(i) = p2_temp(end);
    error_xbis(i) = error_xbis_temp(end);
    error_pbis(i) = error_pbis_temp(end);
    
    error_x(i) = sqrt(x2(i) - (x(i)*x(i)));
    error_p(i) = sqrt(p2(i) - (p(i)*p(i)));
end

%% Figures %%
%%%%%%%%%%%%%

xbis = x;
xbis(end) = [];
dtbis = [0.0625 0.125 0.25 0.5];
[a,b,c,d] = linearFit(dtbis.*dtbis, xbis);
z = polyval([a, b], [0, 1]);
extend = polyval([a, b], [0, 1]);
s = sprintf('Linear fit: y = %.3fx + %.3f',a,b);
figure
plot([0, 1], extend, 'linewidth', lw)
hold on
plot(dt.*dt, x, 'ko', 'MarkerSize', ms)
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('<x>(t_{fin})')
legend(s, 'Data obtained', 'Location', 'nw')
grid on

error_x2 = error_xbis;
error_x2(end) = [];
dtbis = [0.0625 0.125 0.25 0.5];
[a,b,c,d] = linearFit(dtbis.*dtbis, error_x2);
z = polyval([a, b], [0, 1]);
extend = polyval([a, b], [0, 1]);
s = sprintf('Linear fit: y = %.3fx + %.3f',a,b);
figure
plot([0, 1], extend, 'linewidth', lw)
hold on
plot(dt.*dt, error_xbis, 'ko', 'MarkerSize', ms)
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('<\Delta x>(t_{fin})')
legend(s, 'Data obtained', 'Location', 'sw')
grid on

pbis = p;
pbis(end) = [];
dtbis = [0.0625 0.125 0.25 0.5];
[a,b,c,d] = linearFit(dtbis.*dtbis, pbis);
z = polyval([a, b], [0, 1]);
extend = polyval([a, b], [0, 1]);
s = sprintf('Linear fit: y = %.3fx + %.3f',a,b);
figure
plot([0, 1], extend, 'linewidth', lw)
hold on
plot(dt.*dt, p, 'ko', 'MarkerSize', ms)
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('<p>(t_{fin})')
legend(s, 'Data obtained', 'Location', 'se')
grid on

error_p2 = error_pbis;
error_p2(end) = [];
dtbis = [0.0625 0.125 0.25 0.5];
[a,b,c,d] = linearFit(dtbis.*dtbis, error_p2);
z = polyval([a, b], [0, 1]);
extend = polyval([a, b], [0, 1]);
s = sprintf('Linear fit: y = %.5fx + %.5f',a,b);
figure
plot([0, 1], extend, 'linewidth', lw)
hold on
plot(dt.*dt, error_pbis, 'ko', 'MarkerSize', ms)
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('<\Delta p>(t_{fin})')
legend(s, 'Data obtained', 'Location', 'nw')
grid on
