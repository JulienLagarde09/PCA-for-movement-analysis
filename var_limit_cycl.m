 function [data_perp_local_1,data_perp_local_2,var_alongCyc_PCA_1,var_alongCyc_PCA_2,pos_vites,eigvect, eigval,nb_sect] = var_limit_cycl(position,vitesse,wind)

code to estimate variability of noisy self sustained oscillators
along the limit cycle (tangent) and normal to the limit cycle

nb : le nb de point pour estimer var traj normales au cycle = cst


observed variability is considered as the product of the interplay of
stability (determinist) and noise (random, stochastic)
much as defined by drift and diffusion part in a Langevin and Fokker
Planck formulation (see Gardiner, 1993, Haken, 1977, Schoner et al 1986)
hypothesis of additive noise (Gaussian increments, delta correlated)
check if multiplicative noise changes the analysis

polar coord pour estimation variance along direction,
normale au cycle limite
prend N valeurs successives selon ordre de grandeur croissant
de l'angle (theta), en utilisant (SORT) tous cycles confondus,
tous essais confondus, on se retrouve avec les valeurs des angles
prôches ordonnées, ex : 0, 0.1, 0.11, 0.012 etc, 
et prend les valeurs de l'amplitude (rho) correspondantes
et calc la variance et la moyenne

ceci correspond à une découpe de l'histogramme 2D;
en secteurs angulaires du cadran, la longueur du secteur
est réglée par la taille de la fenêtre;
on prend les points contenus dans le premier secteur, et calcule
les stats de leur amplitude

Les isochrones (lignes de phases équivalentes autour
du cycle limite) ne sont pas des lignes droites pour des oscillations
de relaxation (voir par exemple définition de la phase par
Kuramoto (1984)

% INPUTS
wind = wind;
% wind : nb de valeurs dans chaque section
% position
% vitesse

% OUTPUTS
% variance normale au cycle limite : var_alongCyc_PCA
% nombre de sections/cadrans : nb_sect
% données projetées sur le premier composant PCA: data_perp_local
% matrice position et vitesse par section : pos_vites
% output de la PCA sur les points [positions vitesse] par section: eigvect,
% eigval
position = position';
vitesse = vitesse';
[THETA,RHO] = cart2pol(position,vitesse);
a = [THETA; RHO]; % coord polaires du cycle limite
a = a';
a(:,3) =  position;
a(:,4) =  vitesse;


% 3D histogram
% % figure
% % hist3(dat,[100 100]) % Draw histogram for 2D data

% première colonne: theta, autres colonnes les valeurs de rho
% correspondantes
v = sortrows(a,1);% class croissant par theta prem col, donc - pi à pi
vv = v(:,2);% toutes val rho
L = length(vv);% toutes val

% Si prend que valeurs positives (half cycle):
b = find(v(:,1)>0);
vv2 = vv(b(1):end);% val >0
L2 = length(vv2);% val >0

% Si prend que valeurs negative (half cycle):
c = find(v(:,1)<0);
vv3 = vv(1:c(end));% val >0
L3 = length(vv3);% val >0

 nb_sect = floor(L/wind); % nb de sections/cadrans de (n = wind) values
 
% method using PCA

for k = 1:nb_sect
    pos_vites(:,k,:) = v(wind*k - wind +1 : wind*k,3:4);% prend les colonnes 3 et 4, soit x et dx/dt    
end

% CHECK QUADRANTS
% figure
% for i = 1:85
%     plot(pos_vites(:,i,1),pos_vites(:,i,2),'.k')
%     hold on
% end   


for k = 1:nb_sect
[eigvect(:,:,k), scores(:,:,k), eigval(:,k)] = pca(squeeze(pos_vites(:,k,:)));
end
% projeter sur V en ne gardant qu'une direction propre
for k = 1:nb_sect
 data_perp_local_1(:,k) = eigvect(1,:,k)*(squeeze(pos_vites(:,k,:)))';
end
for k = 1:nb_sect
 data_perp_local_2(:,k) = eigvect(2,:,k)*(squeeze(pos_vites(:,k,:)))';
end

var_alongCyc_PCA_1 = var(data_perp_local_1,1);% variance
var_alongCyc_PCA_2 = var(data_perp_local_2,1);% variance
% mean_alongCyc_PCA = mean(data_perp_local,1);% mean


% figure(1)
% plot(var_alongCyc)
% title('variance amplitude(rho)ds windows 500 pts, sorted by theta') 
% % % % % figure(2)
% % % % % hist3(matrix_pos_Vit,[50 50]);% histo 2D
% % % % % view(-17,72);
% % % % % figure(3)
% % % % % bin = 0:1:20;
% % % % % hist(vv,bin)
% % % % % title('histo amplitudes polar raw')
% figure(4)
% plot(mean_alongCyc)
% title('mean amplitude wind')
% % 
% % figure(2)
% % plot(movingvar(v(:,2),200))
% % 
% % figure(3)
% % plot(movingvar(vv,200))
% % % % % figure(6)
% % % % % xx = 1:length(var_alongCyc);
% % % % % [AX,H1,H2] = plotyy(xx,var_alongCyc,xx,mean_alongCyc,'plot');
% % % % % set(get(AX(1),'Ylabel'),'String','var') 
% % % % % set(get(AX(2),'Ylabel'),'String','mean')
% % % % % xlabel('Time') 
% % % % % title('var and mean amplitude limit cycle polar raw')
% % % % % set(H1,'LineStyle','-')
% % % % % set(H2,'LineStyle','-')



