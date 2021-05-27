%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output';
nom = 'Evolution.gif';

data = load([fichier,'_observables.out']);
t = data(:,1);

data = load([fichier,'_psi2.out']);
a = size(t);
psi = cell(1, a(1));
for i=1:a
    psi{i} = data(i+1,:);
end
x=data(1,:);
bite = size(x)
size(psi{1})
h = figure;
for i = 1:size(t)
    %Dessin de la fonction
    plot(x, psi{i}.')
    xlim([-200, 200])
    ylim([-0.2 0.2])
    xlabel('X', 'fontsize', fs)
    ylabel('Y', 'fontsize', fs)
    drawnow
    
    frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im,256); 

  if i == 1 
      imwrite(imind,cm,nom,'gif', 'Loopcount',inf); 
  else 
      imwrite(imind,cm,nom,'gif','WriteMode','append'); 
  end 

end
