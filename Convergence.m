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

repertoire = '~/MATLAB/Exercice8_2021/'; % Chemin d'acces au code compile
executable = 'Exercice8'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

nsimul = 5; % Nombre de simulations a faire

dt = [16 8 4 2 1];


paramstr = 'dt'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul);

for i = 1:nsimul
    filename = [paramstr, '=', num2str(param(i))];
    output{i} = [filename, '.out'];
    c = sprintf('%s%s %s%s %s=%.15g output_observables=%s', repertoire, executable, repertoire, input, paramstr, param(i), output{i})
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

figure
plot(dt.*dt, x, '-k+')
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('Error on <x>')
grid on

figure
plot(dt.*dt, error_x, '-k+')
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('Error on <\Delta x>')
grid on

figure
plot(dt.*dt, p, '-k+')
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('Error on <p>')
grid on

figure
plot(dt.*dt, error_x, '-k+')
set(gca,'fontsize',fs)
xlabel('\Delta t^2')
ylabel('Error on <\Delta x> bis')
grid on
