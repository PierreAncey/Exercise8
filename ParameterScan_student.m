% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '~/MATLAB/Exercice8_2021/'; % Chemin d'acces a l'executable Sur Linux ou MacOSX
% repertoire = ''; % Pour Windows
executable = 'Exercice8'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree

tfin = 2000; % Doit etre le meme que dans le fichier d'entree
dt = [16 8 4 2 1]; % A modifier selon vos besoins!
dt = tfin./round(tfin./dt); % Pour etre sur d'arriver exactement a tfin
nsimul = length(dt); % Nombre de simulations a faire
paramstr = 'dt'; % Nom du parametre a scanner, par exemple dt, w, x0, etc
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output_observables = cell(1, nsimul);

for i = 1:nsimul
    output_observables{i} = ['output_observables_',paramstr, '_', num2str(param(i)),'.out'];
    eval(sprintf('!%s%s %s %s=%.15g output_observables=%s', repertoire, executable, input, paramstr, param(i), output_observables{i}))
    disp('Done.')
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dt'))
    xfin = zeros(1,nsimul);
    dxfin = zeros(1,nsimul);
    pfin = zeros(1,nsimul);
    dpfin = zeros(1,nsimul);
else
    E0 = zeros(1,nsimul);
    trans = zeros(1,nsimul);
end

for ii = 1:nsimul
    data = load(output_observables{ii});
    if(strcmp(paramstr,'dt'))
        xfin(ii) = data(end,6);
        dxfin(ii) = data(end,11);
        pfin(ii) = data(end,8);
        dpfin(ii) = data(end,12);
    else
        E0(ii) = data(1,4);
        trans(ii) = data(end,3);
    end
end

if(strcmp(paramstr,'dt'))
    figure
    subplot 221
    plot(dt.^2,xfin,'k+')
    grid on
    xlabel('\Deltat^2')
    ylabel('<x> (t_{fin})')
    subplot 222
    plot(dt.^2,dxfin,'k+')
    grid on
    xlabel('\Deltat^2')
    ylabel('<\Deltax> (t_{fin})')
    subplot 223
    plot(dt.^2,pfin,'k+')
    grid on
    xlabel('\Deltat^2')
    ylabel('<p> (t_{fin})')
    subplot 224
    plot(dt.^2,dpfin,'k+')
    grid on
    xlabel('\Deltat^2')
    ylabel('<\Deltap> (t_{fin})')
else
    figure
    plot(n,trans)
    grid on
    xlabel('E / V_0')
    ylabel('P_{x>0}(t_{fin})')
end